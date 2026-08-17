# -*- coding: utf-8 -*-
"""M6 临床计算纯函数：eGFR-Schwartz 与 KDIGO 儿童 CKD 分期。

不依赖 fastmcp，可直接 import 单测。所有公式口径在 ADR-006 锁定。

公式来源与口径（ADR-006）：
- eGFR（床旁 Schwartz 2009）：eGFR = 0.413 × 身高(cm) / Scr(mg/dL)，结果以 ml/min/1.73m2 计。
  这是当前 KDIGO / 中国儿科 CKD 随访最常用的体表面积标准化估算式。
- eGFR（含 BUN 修订 Schwartz 2009）：当提供 BUN 且怀疑 eGFR<60 时采用
  eGFR = 0.413 × 身高 / (Scr + 0.003×BUN − 0.024)，对低 eGFR 儿童更准确。
- eGFR（经典 k 值 Schwartz）：eGFR = k × 身高 / Scr，k 随年龄/成熟度的默认取值：
  早产儿 0.33、足月儿(<1y) 0.45、儿童(1–12y) 0.55、青少年(≥13y) 0.70；允许显式覆盖 k_value。
- CKD 分期（KDIGO 2024 儿童）：仅按 eGFR 分 G1–G5；白蛋白尿按 UACR(mg/g) 或 UPCR(mg/g) 分 A1–A3；
  合并分期写作 GxAx（如 G3aA2）。<2 岁婴儿的 eGFR 阈值与年长儿不同，已在 note 中警示。
"""
from __future__ import annotations

import json
import math
import os
import threading

from types import MappingProxyType
from typing import Any, Dict, List, Literal, Mapping, Optional

from a207_policy import enforce_read, get_caller

MCP_NAME = "CKDNutri-assessment-mcp"

# M-1（2026-08-15，第六轮审查）：revised2009 是**无文献出处的自造线性式**——
# 0.413×Ht/(Scr+0.003×BUN−0.024)。Schwartz 2009 含 BUN 的修订式是 CKiD 组合式
# （幂函数且**必需胱抑素 C**）：39.8×[Ht/Scr]^0.456×[1.8/CysC]^0.418×[30/BUN]^0.079
# ×[1.076]male×[Ht/1.4]^0.179（2012 版；2009 原版 39.1×...^0.516×...^0.294×...^0.169
# ×[1.099]male×...^0.188）。本系统无胱抑素 C 数据契约，无法实现真实 CKiD 式——
# **移除假公式**：提供 BUN 时不再自动默认 revised2009，BUN 不参与计算（有告警）。
EgfrMethod = Literal["bedside2009", "classic"]

_BEDSIDE_K = 0.413  # Schwartz 2009 床旁系数
# 肌酐单位换算单一事实源：µmol/L → mg/dL（P1 get_labs 输出 scr_umol_L 转 P4 公式所需单位）
_SCR_UMOL_TO_MGDL = 1.0 / 88.4

# 临床合理性范围（儿童 CKD 系统，超物理区间直接拒绝，防荒谬 eGFR 被分期）
_MAX_AGE_YEARS = 18.0      # 本系统为儿童 CKD（PRNT/KDIGO 儿童适用域）
_MAX_HEIGHT_CM = 250.0     # 儿童身高物理上限（超此值必为录入错误）
_MIN_SCR_MGDL = 0.05       # 血肌酐 mg/dL 下限（低于此值 eGFR 荒谬放大，必为录入错误）

# 规则 schema 枚举与必填键（fail-closed 校验单一事实源，集中模块顶部）
_VALID_RULE_LEVELS = frozenset({"L1", "L2", "L3"})
_VALID_ABS_OPS = frozenset({"gt", "gte", "lt", "lte", "between"})
_VALID_TREND_DIRECTIONS = frozenset({"up", "down"})
_RULE_REQUIRED_KEYS = ("id", "name", "level", "type", "metric", "unit", "description")
_LEVEL_RANK = {"L1": 3, "L2": 2, "L3": 1, "none": 0}


def _egfr_to_g(egfr: float, dialysis_mode: Optional[str] = None) -> str:
    """KDIGO 2024 儿童 eGFR 分期（ml/min/1.73m2）。

    F1 裁决（2026-08-17，十二审）：**G5D 按透析状态定义，与 eGFR 数值无关**——
    KDIGO/NICE/PRNT 把透析依赖患儿单独列为 G5D（随访每月一次），透析患儿即使
    残余肾功能尚可（eGFR 15-30 常见）也应按 G5D 管理。此前 dialysis_mode 检查在
    egfr>=15 数值分期之后，透析患儿被归 G4/G3b（随访 60 天 vs G5D 30 天，
    节奏减半）。现透析检查提前到数值分期之前（约 2 行重排）。
    注：S1（2026-08-13）原注释"eGFR 未跌入 G5 前仍按数值分期"是有意设计，但
    与 KDIGO 实践相悖，本轮裁决修正（临床安全方向：透析随访更密）。
    """
    # F1/F3（2026-08-17）：dialysis_mode 白名单校验提前到数值分期之前——① 透析
    # 状态优先判 G5D（无论 eGFR）；② 录入错误（"hemodialysls" 等）在**任何** eGFR
    # 下都显式拒绝（此前仅 eGFR<15 才触发 raise，eGFR≥15 时被静默忽略 = fail-open）。
    # P4-3（2026-08-15）：白名单校验——此前任意非空非 "none" 字符串（"yes"/
    # "NoNe "/"0" 等录入错误）都被判 G5D（fail-open）。server 层已有 Literal 枚举
    # 拦截，但 core 是纯函数可被编排层直调（绕过 server），补同口径白名单：仅显式
    # 透析模式才判 G5D，其余拒绝（录入错误不应静默升为透析分期）。
    if dialysis_mode:
        dm = dialysis_mode.strip().lower()
        if dm == "none":
            pass
        elif dm in ("hemodialysis", "peritoneal"):
            return "G5D"
        else:
            raise ValueError(
                f"dialysis_mode 必须是 none / hemodialysis / peritoneal 之一，收到："
                f"{dialysis_mode!r}")
    if egfr >= 90:
        return "G1"
    if egfr >= 60:
        return "G2"
    if egfr >= 45:
        return "G3a"
    if egfr >= 30:
        return "G3b"
    if egfr >= 15:
        return "G4"
    return "G5"


def _acr_to_a(uacr_mg_g: float) -> str:
    """白蛋白尿 A 期：UACR mg/g。"""
    if uacr_mg_g < 30:
        return "A1"
    if uacr_mg_g <= 300:
        return "A2"
    return "A3"


def _normalize_sex(sex: Optional[str]) -> Optional[str]:
    """统一清洗 sex 参数：去空格/大小写归一 → "female" | "male" | None。

    所有 sex 判定（_schwartz_k / sex_note / DAG 入口）统一走本函数，防双轨漂移——
    某处只 lower 未 strip 会把 ≥13y 女性静默降级 k=0.70（eGFR 高估 27%）。
    未知值返回 None（向后兼容保持 0.70）。
    """
    if sex is None:
        return None
    s = str(sex).strip().upper()
    if s in ("F", "FEMALE", "女"):
        return "female"
    if s in ("M", "MALE", "男"):
        return "male"
    return None


def _schwartz_k(age_years: float, k_value: Optional[float], is_preterm: bool = False,
                sex: Optional[str] = None) -> float:
    """Schwartz 经典 k 值：早产儿(<37周) 0.33、<1yr 足月儿 0.45、1-12yr=0.55、
    ≥13yr 青少年男性 0.70 / 女性 0.55（经典 Schwartz 1976/1984 参数表）。

    k_value 显式传入优先。早产修正与 ≥13y 性别分化是关键临床项（用错可致 eGFR
    高估 27-36% 并误分期），调用方须按病史与性别正确传参。
    """
    if k_value is not None:
        return float(k_value)
    if age_years < 1:
        return 0.33 if is_preterm else 0.45
    if age_years < 13:
        return 0.55
    if _normalize_sex(sex) == "female":
        return 0.55
    return 0.70


def calc_egfr_schwartz(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    method: EgfrMethod = "bedside2009",
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    serum_creatinine_unit: str = "mg_dL",
    is_preterm: bool = False,
    sex: Optional[str] = None,
) -> dict:
    """估算肾小球滤过率（Schwartz 系列）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    返回：egfr(ml/min/1.73m2)、method、formula(所用公式串)、note(警示/口径)。
    BUG-34（2026-08-12）：显式取 caller 并回写审计字段（此前仅 enforce_read 内部取用，
    调用者身份不落返回，无法追溯谁触发了计算）。
    BUG-40（2026-08-12）：新增 serum_creatinine_unit 单位归一化——P1 clinical-data
    `get_labs` 返回 `scr_umol_L`（µmol/L），本函数历史上只接受 mg/dL；编排层若把 P1
    数值直接透传会错 88.4 倍（eGFR 缩小 88 倍）。现支持 unit="umol_L" 自动 ÷88.4 转
    mg/dL（默认仍为 mg_dL，向后兼容）。也提供显式转换函数 `scr_umol_to_mgdl`。
    BUG-47（2026-08-12）：新增 is_preterm——早产儿经典 k=0.33（仅 <1yr 生效），
    防 eGFR 高估 36%。
    2026-08-12（系统性审查，P1）：新增 sex（M/F/male/female/男/女）——经典式 ≥13y
    青少年女性 k=0.55（男性 0.70），防 eGFR 高估 27%；仅 method="classic" 生效。
    """
    caller = get_caller()
    enforce_read(MCP_NAME)
    # 2026-08-12（系统性审查，P3）：sex 入口归一化一次，sex_note 复用归一化值
    # （_schwartz_k 内部保留 _normalize_sex 防御直接调用者，幂等无副作用）。
    sex = _normalize_sex(sex)
    if method not in ("bedside2009", "classic"):
        # M-1（2026-08-15）：revised2009 已移除（无文献出处的自造线性式）
        if method == "revised2009":
            raise ValueError(
                "method='revised2009' 已移除：该线性式无文献出处。Schwartz 2009 含 BUN "
                "修订为 CKiD 组合式（幂函数且必需胱抑素 C），本系统未实现；请使用 "
                "bedside2009（BUN 不参与）或 classic")
        raise ValueError(
            f"无效的 method: '{method}'，必须为 bedside2009 / classic 之一"
        )
    # BUG-58：显式拒绝 NaN/Inf 入参（NaN 身高/肌酐会算出 NaN eGFR，再被 classify_ckd 静默归 G5）
    # 2026-08-13（policy 审查，高优先级#5）：补**临床合理性上限**——此前仅做有限性与
    # 非负校验，age=200 / height=9999 / scr=1e-6 会算出荒谬 eGFR 并被分期。物理区间
    # 超限一律拒绝（录入错误），常量见模块顶部 _MAX_AGE_YEARS/_MAX_HEIGHT_CM/_MIN_SCR_MGDL。
    age_years = _require_finite(age_years, "age_years")
    height_cm = _require_finite(height_cm, "height_cm")
    if age_years < 0:
        raise ValueError("age_years 不能为负")
    if age_years > _MAX_AGE_YEARS:
        # P4-1（2026-08-15）：硬拒绝（儿童系数用于成人会系统性偏估，fail-closed）+
        # 错误信息直接指引成人公式——此前"提示改用 CKD-EPI"写在返回值 caveat 里，
        # 但 >18 已被此处拒绝、caveat 不可达（死代码），提示从未到达用户。
        raise ValueError(
            f"age_years 超出儿童 CKD 适用域（> {_MAX_AGE_YEARS:.0f} 岁）：本系统仅覆盖"
            "儿童 CKD（Schwartz 公式），成人请改用 CKD-EPI 公式评估")
    if height_cm <= 0:
        raise ValueError("height_cm 必须 > 0")
    if height_cm > _MAX_HEIGHT_CM:
        raise ValueError(f"height_cm 超出物理合理范围（> {_MAX_HEIGHT_CM:.0f} cm）")
    scr = _require_finite(_normalize_scr(serum_creatinine_mgdl, serum_creatinine_unit),
                          "serum_creatinine")
    if scr <= 0:
        raise ValueError("serum_creatinine 必须 > 0")
    if scr < _MIN_SCR_MGDL:
        raise ValueError(f"serum_creatinine 低于物理合理下限（< {_MIN_SCR_MGDL} mg/dL）")
    # M-1（2026-08-15）：revised2009 前置校验随假公式一并移除（method 校验已在上方拒绝）
    # 2026-08-12（系统性审查，P1）：k_value 缺校验——k_value=0 会算出 eGFR=0 → classify
    # 误判 G5；负值产生负 eGFR。显式强校验（与 age/height/scr 同口径）。
    if k_value is not None:
        k_value = _require_finite(k_value, "k_value")
        if k_value <= 0:
            raise ValueError("k_value 必须 > 0")

    if method == "classic":
        k = _schwartz_k(age_years, k_value, is_preterm=is_preterm, sex=sex)
        egfr = k * height_cm / scr
        formula = f"eGFR = k×height/Scr, k={k}"
        k_note = "早产儿 k=0.33" if (is_preterm and age_years < 1) else "经典年龄带 k"
        sex_note = ("；≥13y 女性 k=0.55" if (age_years >= 13 and k_value is None
                                            and sex == "female") else "")
        note = f"{k_note}；<1y=0.45, 1-12y=0.55, ≥13y 男0.70/女0.55，可被 k_value 覆盖。{sex_note}"
        # M-1（2026-08-15）：revised2009 已移除——含 BUN 修订需胱抑素 C（CKiD 组合式）
        if bun_mg_dl is not None:
            note += ("（提供了 BUN 但经典式未使用；含 BUN 的修订需胱抑素 C（CKiD 组合式），"
                     "本系统未实现，BUN 不参与计算）")
    else:  # bedside2009
        egfr = _BEDSIDE_K * height_cm / scr
        formula = "eGFR = 0.413×height/Scr"
        note = "床旁 Schwartz 2009（KDIGO 推荐默认式）。"
        # 2026-08-12（系统性审查，P1）：床旁式固定 k=0.413 无早产修正——若调用方显式
        # 传 method="bedside2009" 且 is_preterm（<1y），k=0.413 相对早产 k=0.33 高估
        # eGFR 约 25%（0.413/0.33）。显式提示改用 classic，不静默吞参。
        if is_preterm and age_years < 1:
            note += (" 注意：is_preterm=True 时床旁式仍用固定 k=0.413（无早产修正），"
                     "相对早产 k=0.33 高估 eGFR 约 25%；如需 k=0.33 请用 method='classic'。")
        # M-1（2026-08-15）：revised2009 已移除——含 BUN 修订需胱抑素 C（CKiD 组合式）
        if bun_mg_dl is not None:
            note += ("（提供了 BUN 但床旁式未使用；含 BUN 的修订需胱抑素 C（CKiD 组合式），"
                     "本系统未实现，BUN 不参与计算）")

    # P4-3（2026-08-15）：eGFR 结果物理上限——入参有下限（scr≥0.05）但无结果上限，
    # 边界组合（如身高 250cm + scr 0.05）可算出 eGFR 数千被静默分期 G1。生理上
    # 儿童 eGFR 罕见 >200（正常 90-140），超限必为录入错误，拒绝（fail-closed）。
    # F1（2026-08-16，第七轮审查）：eGFR 上限 200→250——200 会误拒正常变异（矮小+
    # 低肌酐儿童如 3 岁 100cm + scr 0.2 → eGFR≈206，健康儿童高值可达 200+）；250
    # 以上才是必然的录入错误/单位错配（儿童生理上限 ~230）。
    if egfr > 250:
        raise ValueError(
            f"eGFR 计算结果 {egfr:.1f} ml/min/1.73m² 超出儿童生理合理上限（>250），"
            "通常是身高/肌酐录入错误或单位错配，请核查数据")
    egfr_rounded = round(egfr, 1)
    # N1 修复（2026-08-13）：**未 round 值判级、round 值展示**——此前 round(egfr,1)
    # 先于分期，89.96→90.0 会在边界翻 G2→G1（DAG 用同一舍入值喂 R-07 同理）。
    # 现在返回 egfr_raw（原始，供分期/趋势规则判级）+ egfr（round 展示），
    # 边界不再翻转；note 仍显式标注展示舍入口径供审计。
    note += " eGFR 已四舍五入至 0.1 ml/min/1.73m2（展示值；判级使用未舍入原始值）。"
    pediatric_caveat = (
        "注意：<2 岁婴儿 eGFR 参考范围低于年长儿（正常可低至 60–90），"
        "G1/G2 阈值在婴儿期需结合月龄与生长曲线解读。"
        if age_years < 2 else ""
    )
    # P4-1（2026-08-15）：adult_caveat 恒为空——>18 岁已在入参校验硬拒绝
    # （错误信息含"成人请改用 CKD-EPI"指引），旧版"age>18 提示改用 CKD-EPI"
    # 分支在此处不可达（死代码）。字段保留以维持返回契约（客户端兼容），值恒 ""。
    adult_caveat = ""
    # BUG-15：成功返回统一 {ok, data} 信封（此前扁平结构）；BUG-34：含 caller 审计字段
    return {
        "ok": True,
        "data": {
            "caller": caller,
            "egfr": egfr_rounded,       # 展示值（round 0.1）
            "egfr_raw": egfr,           # N1：未舍入原始值（分期/趋势规则判级用）
            "unit": "ml/min/1.73m2",
            "method": method,
            "formula": formula,
            "pediatric_caveat": pediatric_caveat,
            "adult_caveat": adult_caveat,  # 四审：>18 岁超出儿童公式适用域提示
            "note": note,
            "scr_unit_used": _unit_label(serum_creatinine_unit),
        },
    }


def scr_umol_to_mgdl(value_umol_L: float) -> float:
    """显式单位转换：µmol/L → mg/dL（换算系数单一事实源 _SCR_UMOL_TO_MGDL）。

    供编排层把 P1 get_labs 的 `scr_umol_L` 转成 P4 eGFR 公式所需单位。
    入参 fail-closed：NaN/Inf/负值拒绝（物理上不可能，显式报错）。
    """
    v = _require_finite(value_umol_L, "value_umol_L")
    if v < 0:
        raise ValueError("value_umol_L 不能为负")
    return v * _SCR_UMOL_TO_MGDL


def _normalize_unit(unit: str | None) -> str:
    """单位字符串归一化单一事实源：小写 + 去空格/斜杠 + µ→u，缺省 mg_dL。

    2026-08-13（二审 #4）：兼容带斜杠/空格的写法（"umol/L"、"mg/dL"、"umol/l"）——
    docstring 与文档均用 "µmol/L" 写法，core 是可直接 import 的纯函数库，
    编排层照抄注释写法直调会误抛 ValueError；归一化后再匹配。
    """
    u = (unit or "mg_dL").strip().lower().replace("µ", "u")
    return u.replace("/", "").replace(" ", "")


def _normalize_scr(value: float, unit: str) -> float:
    """按声明单位归一化肌酐到 mg/dL；非法单位显式报错（fail-closed）。

    fail-closed 说明：先 _require_finite 再换算——None/NaN/Inf 统一 ValueError
    （此前 float(None) 抛 TypeError 冒泡，core 直调不优雅）。
    """
    v = _require_finite(value, "serum_creatinine")
    u = _normalize_unit(unit)
    if u in ("mg_dl", "mgdl", "mg"):
        return v
    if u in ("umol_l", "umoll", "umol"):
        return v * _SCR_UMOL_TO_MGDL
    raise ValueError(f"无效的 serum_creatinine_unit: {unit!r}，可用 mg_dL / umol_L")


def _unit_label(unit: str) -> str:
    u = _normalize_unit(unit)
    return "umol/L（已 ÷88.4 转 mg/dL）" if u in ("umol_l", "umoll", "umol") else "mg/dL"


def _require_finite(value: Any, name: str) -> float:
    """把输入转 float 并显式拒绝 NaN/Inf（fail-closed）。

    原因：NaN 静默穿透所有比较（nan<0=False 绕过负数校验、_egfr_to_g(nan) 落 G5），
    inf 恒真落 G1——数值异常必须显式报错，不得静默误分期（临床安全）。
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} 必须为有效数值")
    if math.isnan(v) or math.isinf(v):
        raise ValueError(f"{name} 必须为有效的有限数值")
    return v


def classify_ckd(
    egfr: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    dialysis_mode: Optional[str] = None,
    upcr_mg_mmol: Optional[float] = None,
) -> dict:
    """KDIGO 2024 儿童 CKD 合并分期。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    egfr 必需；白蛋白尿二选一（uacr 或 upcr）。返回 g、a、stage(GxAx)、description、risk_note。

    S1 修复（2026-08-13）：新增 dialysis_mode 可选入参——透析患儿 eGFR<15 分期为
    G5D（KDIGO 2024 / PRNT 2025），与 G5 区分（随访频率不同：G5D 每月 vs G5 每 60 天）。
    合法值：none / hemodialysis / peritoneal（透传 classify 后的 G5D 供 care 层消费）。

    M-5（2026-08-16，十一审）：**UPCR 单位防呆**——P1 契约 upcr_mg_mmol（mg/mmol）
    与 KDIGO UPCR 阈值单位（mg/g）差 **8.84 倍**，此前仅收 upcr_mg_g，编排层误传
    P1 值即静默错 8.84 倍。现支持 upcr_mg_mmol（自动 ÷8.84 换算为 mg/g 后参与）：
    - 仅传 upcr_mg_g：按原口径（KDIGO mg/g）；
    - 仅传 upcr_mg_mmol：÷8.84 换算并标注来源；
    - **同时传两单位：拒绝**（单位歧义，防静默错 8.84 倍）。
    换算系数 8.84（肌酐 mmol→mg 摩尔质量 113.12 mg/mmol ÷ 1000 × ... 即
    1 mg/mmol = 8.84 mg/g，见 clinical-data reference.UNIT_ALIASES 同口径）。

    BUG-34：显式取 caller 并回写审计字段。
    """
    caller = get_caller()
    enforce_read(MCP_NAME)
    # BUG-58：显式拒绝 NaN/Inf（此前 nan 绕过 <0 校验、_egfr_to_g 静默落 G5）
    egfr = _require_finite(egfr, "egfr")
    if egfr < 0:
        raise ValueError("egfr 必须 >= 0")
    if dialysis_mode is not None and not isinstance(dialysis_mode, str):
        raise ValueError("dialysis_mode 必须为字符串（none/hemodialysis/peritoneal）")
    # M-5：UPCR 双单位同时传入 → 单位歧义拒绝（防编排层传错静默 8.84 倍）
    if upcr_mg_g is not None and upcr_mg_mmol is not None:
        raise ValueError(
            "UPCR 同时提供了 mg/g 与 mg/mmol 两种单位（差 8.84 倍）——请只传一种，"
            "避免单位歧义导致分期错误")
    if upcr_mg_g is not None:
        upcr_mg_g = _require_finite(upcr_mg_g, "upcr_mg_g")
        if upcr_mg_g < 0:
            raise ValueError("upcr_mg_g 不能为负")
    if upcr_mg_mmol is not None:
        upcr_mg_mmol = _require_finite(upcr_mg_mmol, "upcr_mg_mmol")
        if upcr_mg_mmol < 0:
            raise ValueError("upcr_mg_mmol 不能为负")
    if uacr_mg_g is not None:
        uacr_mg_g = _require_finite(uacr_mg_g, "uacr_mg_g")
        if uacr_mg_g < 0:
            raise ValueError("uacr_mg_g 不能为负")

    g = _egfr_to_g(egfr, dialysis_mode)
    a: Optional[str] = None
    albuminuria_source = ""
    albuminuria_note = None
    upcr_used = upcr_mg_g if upcr_mg_g is not None else (
        (upcr_mg_mmol / 8.84) if upcr_mg_mmol is not None else None)
    if uacr_mg_g is not None:
        a = _acr_to_a(uacr_mg_g)
        albuminuria_source = "UACR"
        if upcr_used is not None:
            # 2026-08-12（系统性审查，P3）：双指标同时传入时透明化——KDIGO 2024 以
            # UACR 为准，但 UPCR 的忽略状态须显式标注供上游审计（此前静默丢弃，
            # 审计无法知晓 UPCR 已提供却未参与分期）。
            albuminuria_note = ("同时提供 UACR 与 UPCR：按 KDIGO 2024 以 UACR 为准，"
                                "UPCR 未用于 A 分期。")
    elif upcr_used is not None:
        # BUG-48（2026-08-12）：KDIGO 2024 白蛋白尿 A 分期**仅基于 UACR**。UPCR 反映尿总蛋白
        # （含球蛋白），与白蛋白不等价——肾病综合征等患儿 UPCR 显著高于 UACR，直接映射
        # 会把 A1 误判为 A2/A3。故仅提供 UPCR 时**不再映射 A 分期**（a=None），
        # 返回明确提示，请补充 UACR 后再做白蛋白尿分期。
        a = None
        albuminuria_source = "UPCR"
        albuminuria_note = ("仅提供 UPCR（尿总蛋白/肌酐比）：KDIGO 2024 白蛋白尿 A 分期仅基于 "
                            "UACR，UPCR 含球蛋白排泄、与白蛋白不等价，未映射 A1/A2/A3。"
                            "请补充 UACR 以完成分期。")

    stage = f"{g}{a}" if a else g

    g_desc = {
        "G1": "eGFR ≥ 90，肾功能正常/亢进（需合并结构性病变才定义 CKD）",
        "G2": "eGFR 60–89，轻度下降",
        "G3a": "eGFR 45–59，轻-中度下降",
        "G3b": "eGFR 30–44，中-重度下降",
        "G4": "eGFR 15–29，重度下降",
        "G5": "eGFR < 15，肾衰竭",
        # S1 修复（2026-08-13）：透析患儿单独列 G5D——随访频率与 G5 不同
        "G5D": "eGFR < 15 且正在透析（G5D）",
    }[g]
    a_desc = None
    if a is not None:
        # BUG-67 后补（2026-08-12）：a_desc 仅当 a 已判定时构建——a=None（仅 UPCR 时）
        # 无需白蛋白尿描述字典，移入分支内代码更清爽（行为不变）。
        a_desc = {
            "A1": "白蛋白尿正常-轻度增高",
            "A2": "白蛋白尿中度增高",
            "A3": "白蛋白尿重度增高",
        }

    risk_note = _risk_note(g, a)

    return {
        "ok": True,
        "data": {
            "caller": caller,  # BUG-34 审计字段
            "stage": stage,
            "g": g,
            "a": a,
            "g_description": g_desc,
            # 💭（2026-08-12 五包审查）：简化冗余条件——a_desc 仅当 a 非 None 时构建
            # （BUG-67 后补），`and a_desc` 恒真，保留 `a is not None` 即可。
            "a_description": a_desc.get(a) if a is not None else None,
            "albuminuria_source": albuminuria_source or None,
            "albuminuria_note": albuminuria_note,
            # 2026-08-12（系统性审查，P3）：双指标同时传入时为 True——机器可读标注
            # UPCR 被忽略（与 albuminuria_note 文本互补，供上游审计程序化判断）。
            "upcr_ignored": (uacr_mg_g is not None and upcr_mg_g is not None),
            "risk_note": risk_note,
        },
    }


def _risk_note(g: str, a: Optional[str]) -> str:
    """按 G/A 给出进展风险与随访强度提示（信息性）。"""
    # S1 修复（2026-08-13）：G5D 与 G5 同属最高风险档（rank=5）
    g_rank = {"G1": 0, "G2": 1, "G3a": 2, "G3b": 3, "G4": 4, "G5": 5, "G5D": 5}[g]
    a_rank = {"A1": 0, "A2": 1, "A3": 2}.get(a or "A1", 0)
    # F2（2026-08-17，十二审）：**G 阶段地板**——KDIGO 2024 风险热图把 G4/G5/G5D
    # 定义为高进展风险（红色区），**与白蛋白尿无关**（G4/G5 无论 A1/A2/A3 都是高
    # 风险）。此前纯加性评分 g_rank+a_rank 阈值 ≥6 高：G4 A1=4、G5 A1=5、G5D A1=5
    # 全落"中等"（3-6 月随访）→ 透析/晚期患儿被建议 3-6 月随访且风险文案降级，
    # 临床分级不足。加地板：G4+ 直接判高（与 KDIGO 热图一致）。
    if g in ("G4", "G5", "G5D"):
        note = "进展风险高：建议缩短随访间隔（如 1–3 个月）并由肾科密切管理。"
    else:
        score = g_rank + a_rank
        if score >= 6:
            note = "进展风险高：建议缩短随访间隔（如 1–3 个月）并由肾科密切管理。"
        elif score >= 3:
            note = "进展风险中等：建议 3–6 个月随访一次。"
        else:
            note = "进展风险低：建议 6–12 个月常规随访。"
    # BUG-58（2026-08-12）：a=None（未提供 UACR / 仅提供 UPCR）时评分按无白蛋白尿计，
    # 可能低估进展风险——显式注明白蛋白尿未计入，提示补充 UACR 复评（透明，不捏造分期）。
    if a is None:
        note += ("（白蛋白尿 A 分期未计入评分：未提供 UACR 或仅提供 UPCR；"
                 "实际进展风险可能更高，请补充 UACR 后复评。）")
    return note

_RULES_PATH = os.path.join(os.path.dirname(__file__), "data", "rules.json")
_RULES: Optional[Dict[str, Any]] = None
_RULES_VIEW: Optional[Mapping[str, Any]] = None
# 懒加载并发锁：FastMCP/多 worker 线程并发首次调用时防重复 I/O 与重复 JSON 解析
# （double-checked locking，GIL 下安全）。
_RULES_LOCK = threading.Lock()


def _load_rules() -> Mapping[str, Any]:
    global _RULES, _RULES_VIEW
    if _RULES is None:
        with _RULES_LOCK:
            if _RULES is None:
                try:
                    with open(_RULES_PATH, "r", encoding="utf-8") as f:
                        _RULES = json.load(f)
                except FileNotFoundError:
                    # 2026-08-12（系统性审查）：异常消息只透出文件名（basename）不泄露服务端
                    # 目录结构（与 clinical-data store 同口径）；完整路径仅留在服务端排障。
                    raise FileNotFoundError(
                        f"风险规则文件缺失：{os.path.basename(_RULES_PATH)}；请确认 data/rules.json 存在"
                    )
                except json.JSONDecodeError as e:
                    raise ValueError(
                        f"风险规则文件 JSON 解析失败：{os.path.basename(_RULES_PATH)}，{e}"
                    )
                # 四审（2026-08-12）：加载期规则结构校验（fail-closed，对齐 content
                # _load_guides/_load_sops 的加载校验）——此前 rules.json 缺键/类型错
                # 在运行时 KeyError（归 INTERNAL_ERROR 但无定位）或**静默误判**：
                # level 未知值 rank=0（漏报危急）、operator 未知值 hit=False（漏报）、
                # trend direction 非法值落 down 分支（方向反判）。配置错误必须在
                # 加载期显式拒绝，不进入评估路径。
                # A-B1 修复（2026-08-14）：校验失败必须**重置缓存再抛**——此前
                # `_RULES = json.load(...)` 在校验前赋值，schema 校验失败后 `_RULES`
                # 非 None 但 `_RULES_VIEW` 仍为 None，模块被投毒：后续所有调用跳过
                # 加载直接返回 None（真实配置错误被包装成 INTERNAL_ERROR 且进程内
                # 永不恢复）。现在校验/归一化失败时 `_RULES = None`，下次调用重新
                # 读取并再次校验（修复配置后自动恢复）。
                try:
                    _validate_rules_schema(_RULES)
                except Exception:
                    _RULES = None
                    raise
                # 2026-08-12（系统性审查，P3）：加载时一次性归一化规则阈值类型——
                # absolute/between 的 low/high 统一 float（int 配置 {"low":1} 在此转 1.0），
                # 保证返回引用的展示契约（"[1.0, 2.0)"）在各配置来源下绝对一致，命中路径
                # 无需再操心类型（_eval_rule 保留的 float() 为幂等防御，float(5.0) 返回原对象）。
                for _r in _RULES.get("rules", []):
                    if _r.get("type") == "absolute" and _r.get("operator") == "between":
                        _r["low"] = float(_r["low"])
                        _r["high"] = float(_r["high"])
                # 2026-08-12（系统性审查，P2）：返回 MappingProxyType 顶层只读视图——
                # 替代裸 dict 引用：类型层面表达只读契约，外部对**顶层键**的增删改
                # （_RULES["rules"]=[...] 等）即刻 TypeError。注意嵌套结构（list/dict）
                # 仍可变（Python 无深层不可变），深层防护依赖"公开 API 不暴露引用 +
                # 调用方自律"（evaluate/list_rules/assess 返回值均为新建对象，已满足）。
                _RULES_VIEW = MappingProxyType(_RULES)
    return _RULES_VIEW


def _reset_rules_cache() -> None:
    """测试辅助：重置规则库懒加载缓存（force_reload 语义）。生产路径勿调。"""
    global _RULES, _RULES_VIEW
    with _RULES_LOCK:
        _RULES = None
        _RULES_VIEW = None


def _validate_rules_schema(rules_doc: Dict[str, Any]) -> None:
    """规则库结构校验（fail-closed，四审 2026-08-12）。

    在加载期一次性校验，杜绝运行时三种静默误判：
    - level 非法值 → _LEVEL_RANK.get(...,0) 归 0（危急规则被静默降级为不触发）；
    - absolute operator 非法值 → _eval_rule hit=False（漏报）；
    - trend direction 非法值 → else 分支当 down 处理（方向反判）。
    """
    rules = rules_doc.get("rules")
    if not isinstance(rules, list) or not rules:
        # #8（2026-08-15）：规则库损坏是**服务端数据问题**，抛 RuntimeError 归
        # INTERNAL_ERROR（detail 脱敏）——此前抛 ValueError 被 translate_error
        # 归 INVALID_INPUT（客户端错误码）+ detail 泄露规则全文，与 P5 数据错误
        # 处理方向相反（服务端数据损坏 ≠ 客户端入参错误）。
        raise RuntimeError("风险规则库缺少非空 'rules' 列表，拒绝加载（服务端数据损坏）")
    # P1-6 修复（2026-08-13）：metric 值域校验——合法短名 = _LAB_ALIAS_TO_RULE 的
    # short 集合（scr/k/p/hb/ca/na/ua/egfr/bun/uacr/upcr）。此前拼错的 metric（如
    # "ks"）静默跳过该规则 → overall_level="none" 给临床"无风险"假象（fail-open）。
    _VALID_RULE_METRICS = frozenset(short for _, (short, _f) in
                                    ((k, v) for k, v in _LAB_ALIAS_TO_RULE.items()))
    for r in rules:
        if not isinstance(r, dict):
            raise RuntimeError(f"规则条目 {r!r} 非字典，拒绝加载")
        missing = [k for k in _RULE_REQUIRED_KEYS if not str(r.get(k) or "").strip()]
        if missing:
            raise RuntimeError(
                f"规则 {r.get('id', '?')} 缺少必填键 {missing}，拒绝加载（fail-closed）")
        if r["level"] not in _VALID_RULE_LEVELS:
            raise RuntimeError(
                f"规则 {r['id']} level={r['level']!r} 非法，必须是 {sorted(_VALID_RULE_LEVELS)}")
        # P1-6：metric 必须落在合法短名集合内（fail-closed，防静默跳过规则）
        if r.get("metric") not in _VALID_RULE_METRICS:
            raise RuntimeError(
                f"规则 {r['id']} metric={r.get('metric')!r} 非法，必须是 "
                f"{sorted(_VALID_RULE_METRICS)}（P1-6：拼错 metric 会静默漏报）")
        if r["type"] == "absolute":
            op = r.get("operator")
            if op not in _VALID_ABS_OPS:
                raise RuntimeError(
                    f"规则 {r['id']} operator={op!r} 非法，必须是 {sorted(_VALID_ABS_OPS)}")
            if op == "between":
                low, high = r.get("low"), r.get("high")
                if not (isinstance(low, (int, float)) and isinstance(high, (int, float))
                        and low < high):
                    raise RuntimeError(
                        f"规则 {r['id']} between 需 low < high 数值，收到 low={low!r} high={high!r}")
            elif not isinstance(r.get("threshold"), (int, float)):
                raise RuntimeError(
                    f"规则 {r['id']} 单边 operator 需数值 threshold，收到 {r.get('threshold')!r}")
        elif r["type"] == "trend_pct":
            direction = r.get("direction")
            if direction not in _VALID_TREND_DIRECTIONS:
                raise RuntimeError(
                    f"规则 {r['id']} direction={direction!r} 非法，必须是 up / down")
            # 2026-08-13（二审 #3）：单阈值（threshold_pct）与区间（low_pct/high_pct）
            # **二选一互斥**——此前 R-08 曾同时携带 threshold_pct:30 与 low/high:30/50
            # （矛盾配置被静默接受，语义由"区间优先"掩盖）；fail-closed 拒绝混用。
            has_single = "threshold_pct" in r
            has_range = "low_pct" in r or "high_pct" in r
            if has_single and has_range:
                raise RuntimeError(
                    f"规则 {r['id']} threshold_pct 与 low_pct/high_pct 必须二选一（互斥），"
                    f"当前同时提供")
            if has_range:
                low_pct, high_pct = r.get("low_pct"), r.get("high_pct")
                if not (isinstance(low_pct, (int, float)) and isinstance(high_pct, (int, float))
                        and low_pct < high_pct):
                    raise RuntimeError(
                        f"规则 {r['id']} 区间型趋势需 low_pct < high_pct，收到 "
                        f"low_pct={low_pct!r} high_pct={high_pct!r}")
            elif not isinstance(r.get("threshold_pct"), (int, float)):
                raise RuntimeError(
                    f"规则 {r['id']} 单阈值趋势需数值 threshold_pct，收到 {r.get('threshold_pct')!r}")
        else:
            raise RuntimeError(
                f"规则 {r['id']} type={r['type']!r} 非法，必须是 absolute / trend_pct")


def _pct_change(new: float, old: float) -> float:
    """相对变化百分比（正=升，负=降）。old<=1e-9 或任一值为 NaN/Inf 视为无效。

    BUG-58（2026-08-12）：原 `old is None` 判断位于 float(old) 之后属死代码
    （None 会在 float() 抛 TypeError 提前返回）；且 NaN 穿透 `old <= 0` 比较后
    会带着 NaN 继续计算（虽被上层 pct != pct 过滤），显式拒绝更清晰。
    2026-08-12（系统性审查，P1）：极小正基线（如 o=1e-321 的微量化验值/浮点误差）
    使 (n-o)/o*100 产生 inf 或 1e300 级巨值——上层 `pct != pct` 只拦 NaN、拦不住
    inf，命中趋势规则后 `round(abs(pct),1)` 抛 OverflowError 击穿 fail-closed。
    ① 基线阈值收紧至 o<=1e-9（正常化验值量级远大于此，不受影响）；
    ② 计算结果二次有限性校验（防分母非极小但分子极端导致的 inf/nan）。
    """
    try:
        n, o = float(new), float(old)
    except (ValueError, TypeError):
        return float("nan")
    # 2026-08-12（系统性审查，P1 二修）：基线阈值 1e-9 → 1e-6——o=1e-8（>1e-9）时
    # (n-o)/o*100 仍可达 1e10 级夸张百分比（非 inf，能穿过上轮 isinf 拦截），命中
    # 趋势规则产生 AKI 类误报。o <= 1e-6 同时覆盖**负数基线**（化验值物理上不可能为负，
    # 视为脏值拒绝，防 down 规则误报）；正常化验值量级远大于 1e-6，零误伤。
    if math.isnan(n) or math.isnan(o) or math.isinf(n) or math.isinf(o) or o <= 1e-6:
        return float("nan")
    res = (n - o) / o * 100.0
    if math.isnan(res) or math.isinf(res):
        return float("nan")
    return res


def _eval_rule(rule: Dict[str, Any], new_labs: Dict[str, float],
               prior_labs: Optional[Dict[str, float]]) -> Optional[Dict[str, Any]]:
    """评估单条规则；命中返回观测细节，未命中/无法评估返回 None。"""
    metric = rule["metric"]
    rtype = rule["type"]
    if rtype == "absolute":
        if metric not in new_labs or new_labs[metric] is None:
            return None
        try:
            v = float(new_labs[metric])  # 防御 JSON 序列化字符串类型
        except (ValueError, TypeError):
            return None
        op = rule["operator"]
        if op == "gt":
            hit = v > rule["threshold"]
        elif op == "gte":
            hit = v >= rule["threshold"]
        elif op == "lt":
            hit = v < rule["threshold"]
        elif op == "lte":
            hit = v <= rule["threshold"]
        elif op == "between":
            low, high = rule["low"], rule["high"]
            hit = low <= v < high
        else:
            # schema 校验（_validate_rules_schema）已 fail-closed 拒绝非法 operator，
            # 运行期不可达；保留显式报错（而非静默 hit=False 漏报）以防 schema 漏判。
            raise ValueError(f"规则 {rule.get('id', '?')} 的 operator 非法: {op!r}")
        if not hit:
            return None
        op_label = {"gt": ">", "gte": ">=", "lt": "<", "lte": "<="}.get(op, op)
        return {
            "metric": metric,
            # 2026-08-12（系统性审查）：observed 统一 round 2 位——umol/L→mg/dL 换算
            # 可能产生长浮点（100/88.4=1.131221719...），直出审计/医生端展示体验差。
            "observed": round(v, 2),
            # 2026-08-12（系统性审查）：threshold 统一为 str——此前 between 返回
            # Tuple[low, high]（(1.0,2.0)）而单边/趋势返回 str（">= 1.5"/"up [30,50)%"），
            # 下游强类型解析器（TS/Pydantic DTO）反序列化会因类型漂移崩溃。
            # 2026-08-12（系统性审查，P2）：float() 幂等防御——_load_rules 加载时已把
            # between 的 low/high 归一化为 float，此处 float() 对 float 输入返回原对象
            # （零开销），仅防本函数被直接调用时收到手工构造的 int 配置。
            "threshold": (f"[{float(rule['low'])}, {float(rule['high'])})" if op == "between"
                          else f"{op_label} {rule['threshold']}"),
            "unit": rule["unit"],
        }
    if rtype == "trend_pct":
        if metric not in new_labs or new_labs[metric] is None:
            return None
        if prior_labs is None or metric not in prior_labs or prior_labs[metric] is None:
            # 无历史对照：趋势类规则无法评估，直接跳过（绝不兜底）
            return None
        pct = _pct_change(new_labs[metric], prior_labs[metric])
        if pct != pct:  # NaN
            return None
        direction = rule["direction"]
        if direction == "up":
            if "low_pct" in rule:  # 区间型（如 R-08 30-50%）
                hit = rule["low_pct"] <= pct < rule["high_pct"]
            else:
                hit = pct >= rule["threshold_pct"]
        else:  # down：同样支持区间型（如"下降 30%–50%"）
            if "low_pct" in rule:
                # 下降区间镜像上升 [low, high)：high 端转负后开放，low 端转负后闭合。
                # 例 [30, 50)→(-50, -30]，即下降 30% 命中边界而 50% 不命中。
                # 2026-08-12（系统性审查）：边界契约确认——本实现与 [low,high) 半开区间
                # 语义严格一致：pct=-50（下降 50%）归入下一阶梯（(-70,-50]），本阶梯不命中；
                # pct=-30 恰好命中本阶梯下界。**配置方须全库统一 [low, high) 半开区间**
                # 语义（low 闭合、high 开放）定义相邻下降区间（如 [30,50) 与 [50,70)），
                # 相邻阶梯在临界点无缝衔接（下降 pct=-high 恒归下一阶梯），不会出现
                # "两侧都不命中（漏报）或两侧都命中（误报）"。
                hit = -rule["high_pct"] < pct <= -rule["low_pct"]
            else:
                hit = pct <= -rule["threshold_pct"]
        if not hit:
            return None
        # 生成阈值描述：区间型用 low_pct/high_pct，单阈值型用 threshold_pct
        if "low_pct" in rule:
            thresh_str = (f"{direction} [{rule['low_pct']}%, {rule['high_pct']}%)")
        else:
            thresh_str = (f">= {rule['threshold_pct']}%" if direction == "up"
                          else f"<= -{rule['threshold_pct']}%")
        return {
            "metric": metric,
            # 2026-08-12（系统性审查，P3）：down 方向 observed 取绝对值——_pct_change
            # 下降为负（-35.0），与阈值描述（"down [30%, 50%)"）同框展示时负数直观上
            # "不在区间内"引发 UI/临床困惑。方向语义已由 threshold 文本承载，取绝对值
            # 不丢信息；up 方向 pct 恒正，abs 无影响。
            "observed": round(abs(pct), 1),
            "threshold": thresh_str,
            "unit": "%",
        }
    return None


# --- BUG-45 修复（2026-08-12）-----------------------------------------------
# P1 clinical-data get_labs 返回**完整键名 + PCP 单位**（scr_umol_L µmol/L、k_mmol_L 等），
# rules.json 的 metric 用**短名 + 规则单位**（scr=mg/dL、k=mmol/L、hb=g/L、ua=umol/L…）。
# 此前引擎"不做单位换算、调用方须传规范单位"——编排层把 P1 结果直传会因键名不匹配
# 使全部 absolute 规则静默跳过（overall_level="none" 假象），且 scr 单位差 88.4 倍。
# 现入口统一归一化：完整键名 → 短名 + 单位换算（仅 scr_umol_L 需要 ÷88.4，其余单位一致）。
_LAB_ALIAS_TO_RULE: dict[str, tuple[str, float]] = {
    "scr_umol_L": ("scr", _SCR_UMOL_TO_MGDL),   # µmol/L → mg/dL（单一事实源）
    # 2026-08-13（二审 #2）：补 scr mg/dL 完整键名变体——与 bun_mg_dL/bun_mg_dl 对称
    # （P1 输入模型/别名表存在 scr_mg_dL 大写 L 风格）。此前缺失：编排层直传
    # new_labs={"scr_mg_dL":...} 时 R-01（L1 AKI 危急规则）静默不触发，与 BUG-45
    # 修复初衷相悖。两变体系数 1.0（mg/dL 与规则单位一致）。
    "scr_mg_dl": ("scr", 1.0),
    "scr_mg_dL": ("scr", 1.0),
    "k_mmol_L": ("k", 1.0),
    "p_mmol_L": ("p", 1.0),
    "hb_g_L": ("hb", 1.0),
    "ca_mmol_L": ("ca", 1.0),
    "na_mmol_L": ("na", 1.0),
    "ua_umol_L": ("ua", 1.0),
    "egfr_ml_min": ("egfr", 1.0),
    # BUG-58（2026-08-12）：bun 短名口径与 mg/dL 一致（M-1 后 revised2009 已移除，保留换算防未来 bun 规则偏差）——
    # P1 的 bun_mmol_L（mmol/L 尿素氮）须 ×2.8 换算，原系数 1.0 会在未来新增 bun
    # 规则时产生 2.8 倍偏差（当前规则库尚无 bun 规则，属隐患修复）。
    "bun_mmol_L": ("bun", 2.8),
    # 2026-08-12（系统性审查，P1）：补 bun mg/dL 完整键名变体——P1 输入模型/别名表
    # 存在 bun_mg_dL（reference.py:219），编排层可能把该键直传 new_labs；此前无法映射
    # 为规则短名 "bun" → assess_clinical_status 的 `"bun" in norm_inputs` 补齐逻辑
    # 失效 → method 自动推理看不到 BUN（M-1 后 BUN 本就不参与，仍保留映射供未来规则）。
    # 两变体（全小写 / P1 大写 L 风格）系数均为 1.0（mg/dL 与规则单位一致）。
    "bun_mg_dl": ("bun", 1.0),
    "bun_mg_dL": ("bun", 1.0),
    # 2026-08-12（系统性审查）：补 uacr/upcr 完整键名映射——P1 返回 uacr_mg_g/upcr_mg_g。
    # 当前规则库尚无蛋白尿规则（白蛋白尿 A 分期走 classify_ckd），映射为防未来规则
    # 扩展 + 键名规范化一致性（与 scr/k 等同模式，系数 1.0）。
    "uacr_mg_g": ("uacr", 1.0),
    "upcr_mg_g": ("upcr", 1.0),
}


def _normalize_labs(labs: Optional[Dict[str, float]]) -> Optional[Dict[str, float]]:
    """把 P1 完整键名 + PCP 单位归一化为规则短名 + 规则单位（BUG-45）。

    已用短名的输入原样保留；完整键名在短名缺失时转换并换算单位。
    返回新 dict，不修改调用方对象。

    2026-08-12（系统性审查，P0）：对全部输入值强校验 `_require_finite`——此前 NaN/Inf
    （如 {"k": float("nan")}）会静默穿透比较（nan > 5.5 恒 False）→ 规则漏检 →
    overall_level="none" 给临床"无风险"假象（fail-open，违背 BUG-58 的 fail-closed 原则）。
    非法值抛 ValueError（冒泡 server._invalid 归 INVALID_INPUT），绝不静默放行。
    2026-08-12（系统性审查）：入口类型强校验——LLM/编排层可能把嵌套对象序列化为 JSON
    字符串直传（如 '{"scr": 1.2}'）；FastMCP/pydantic 参数校验层会拦截该形态，但编排层
    **绕过 FastMCP 直接 import core 调用**时（core 为纯函数库），str 会走 labs.items()
    抛 AttributeError（内部 bug 误报）。显式拒绝非 dict（fail-closed ValueError）。
    """
    if labs is None:
        return None
    if not isinstance(labs, dict):
        raise ValueError(
            f"labs 必须为 dict（键→数值映射），实际为 {type(labs).__name__}；"
            "不接受 JSON 字符串等其它形态")
    if not labs:
        return labs
    out: Dict[str, float] = {}
    for k, v in labs.items():
        fv = _require_finite(v, f"labs[{k!r}]")
        # A-B2 修复（2026-08-14）：规则路径物理区间校验（fail-closed，对齐 eGFR 计算
        # 路径）——此前 {"hb": -10} 触发"重度贫血 R-04"、负 eGFR 参与 AKI 判定，
        # 负值/荒谬值物理不可能却静默参与临床判定（同包两套校验深度）。生化浓度
        # 指标一律 ≥ 0；负值显式拒绝（冒泡 server._invalid 归 INVALID_INPUT）。
        if fv < 0:
            raise ValueError(f"labs[{k!r}] 必须 >= 0（生化指标物理上不可能为负），收到 {fv!r}")
        out[k] = fv
    # 2026-08-12（系统性审查，P2）：完整键名匹配大小写容错——`_normalize_scr` 已对单位
    # lower 化，但此处键名此前严格区分大小写：上游/P1 若传 Scr_umol_L、bun_mg_DL 等变体
    # 无法识别 → scr/bun 相关规则静默跳过 → overall_level="none" 假无风险。现用 lower_map
    # 归一化查找（短名仍精确匹配原样保留，不受影响）。
    lower_map = {k.lower(): k for k in out}
    for full_key, (short, factor) in _LAB_ALIAS_TO_RULE.items():
        matched_key = lower_map.get(full_key.lower())
        if matched_key is not None and short not in out:
            # 完整键名已由 _require_finite 校验为有限数，换算必然成功（原 try/except
            # 静默跳过改为强校验后不再需要——脏值在上一循环已被拦截）
            out[short] = round(out[matched_key] * factor, 4)
    # A-B5 修复（2026-08-14）：**短名键也大小写容错**——docstring 声称"键名匹配不区分
    # 大小写"，但此前只对完整键名做 lower 映射，短名（scr/k/hb）大写传入（如 "HB"）→
    # 规则静默跳过 → overall_level="none" 假象。合法短名集合 lower 后归一化到规范短名。
    _short_set = {short for _, (short, _f) in _LAB_ALIAS_TO_RULE.items()}
    for key in list(out.keys()):
        lk = key.lower()
        if lk in _short_set and lk != key:
            out[lk] = out.pop(key)
    return out


def evaluate_risk_rules(
    new_labs: Dict[str, float],
    prior_labs: Optional[Dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> Dict[str, Any]:
    """基于本轮新数据重评风险。

    入参：
      new_labs: 本轮新化验/测量。**接受两种键名**（BUG-45 归一化）：
        - 规则短名：scr(mg/dL), k, p, hb(g/L), ca, ua(umol/L), egfr(ml/min/1.73m^2)
        - P1 完整键名：scr_umol_L(µmol/L), k_mmol_L, p_mmol_L, hb_g_L, ca_mmol_L,
          ua_umol_L, egfr_ml_min, na_mmol_L, bun_mmol_L —— 自动映射到短名并换算单位
        （scr_umol_L 自动 ÷88.4 转 mg/dL；其余单位两域一致直接透传）。
      prior_labs: 上轮同指标值（用于趋势类规则：R-01/R-08 肌酐、R-07 eGFR）。缺则趋势规则不触发。
      prior_level: 历史系统等级（L1/L2/L3）。**仅作对比输出，绝不兜底**。
      身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。

    返回：
      matched_rules: 命中规则明细列表

      overall_level: 最高等级（L1>L2>L3），无命中为 "none"
      prior_comparison: {prior_level, current_level, delta_note}
      level_correction_applied: True（声明范式已执行）
    """
    caller = get_caller()  # BUG-34：显式取 caller 回写审计字段
    enforce_read(MCP_NAME)
    # BUG-45：键名 + 单位归一化（P1 完整键名/µm→mg/dL），防编排层直传静默失效
    new_labs = _normalize_labs(new_labs) or {}
    prior_labs = _normalize_labs(prior_labs)
    rules_doc = _load_rules()
    matched: List[Dict[str, Any]] = []
    for rule in rules_doc["rules"]:
        detail = _eval_rule(rule, new_labs, prior_labs)
        if detail is None:
            continue
        matched.append({
            "id": rule["id"],
            "name": rule["name"],
            "level": rule["level"],
            "description": rule["description"],
            "observed": detail["observed"],
            "threshold": detail["threshold"],
            "unit": detail["unit"],
        })

    # 最高等级（防御性读取：rules.json 中未知等级默认 rank=0，不抛 KeyError）
    overall = "none"
    for m in matched:
        if _LEVEL_RANK.get(m["level"], 0) > _LEVEL_RANK.get(overall, 0):
            overall = m["level"]

    # 2026-08-12（系统性审查，P3）：命中列表按风险等级降序（L1 危急优先）——rules.json
    # 文件顺序非等级序（如 R-01 L1 后紧跟 R-08 L3），医生/前端审阅 matched_rules 时
    # 应先见危急项再见低危项。overall_level 计算与顺序无关，此排序纯展示层无副作用。
    matched.sort(key=lambda m: _LEVEL_RANK.get(m["level"], 0), reverse=True)

    # 四审（2026-08-12）：空输入显式提示——编排层直接调 evaluate_risk_rules（绕过
    # DAG）传空/缺失 new_labs 时，overall_level="none" 会给"已全面评估且无风险"的
    # 假象；与 DAG 的 risk_completeness 同精神，函数级透明化（DAG 场景 new_labs
    # 非空，此分支不触发）。
    evaluation_note = None
    if not new_labs:
        evaluation_note = (
            "未提供任何化验指标（new_labs 为空），规则引擎未评估任何规则；"
            "overall_level=none 不代表无风险，请补充化验数据后重评。")

    # 对比：仅展示，不兜底
    delta_note = "无历史等级对照"
    if prior_level:
        # 2026-08-12（系统性审查）：prior_level 非法值（如 "UNKNOWN"）此前会落
        # _LEVEL_RANK.get(...,0)=0 并输出"较历史等级 UNKNOWN 降至 none"的语义混乱文案——
        # 强制校验合法枚举 {"L1","L2","L3","none"}，非法值忽略对比（不误导）。
        if prior_level not in _LEVEL_RANK:
            delta_note = f"历史等级 {prior_level!r} 无效，已忽略对比"
        elif prior_level == overall:
            delta_note = f"与历史等级 {prior_level} 持平"
        elif _LEVEL_RANK[overall] > _LEVEL_RANK.get(prior_level, 0):
            delta_note = f"较历史等级 {prior_level} 升高至 {overall}"
        else:
            delta_note = f"较历史等级 {prior_level} 降至 {overall}"

    return {
        "ok": True,
        "data": {
            "caller": caller,  # BUG-34 审计字段
            "matched_rules": matched,
            "overall_level": overall,
            "evaluation_note": evaluation_note,  # 四审：空输入/未评估的显式提示
            "prior_comparison": {
                "prior_level": prior_level,
                "current_level": overall,
                "delta_note": delta_note,
            },
            "level_correction_applied": True,
        },
    }


def list_rules() -> Dict[str, Any]:
    """返回规则清单（不含评估逻辑）。

    2026-08-12（系统性审查）：签名修正为 Dict 信封——此前声明 List[Dict] 但实际返回
    {"ok": true, "data": {"rules": [...]}}，类型检查器报错且上游按 List 索引会 KeyError: 0。
    """
    caller = get_caller()
    enforce_read(MCP_NAME)
    rules_doc = _load_rules()
    out = []
    for r in rules_doc["rules"]:
        entry = {"id": r["id"], "name": r["name"], "level": r["level"],
                 "type": r["type"], "description": r["description"], "unit": r["unit"]}
        if r["type"] == "absolute":
            entry["criterion"] = (f"{r['operator']} {r.get('threshold')}"
                                  if r["operator"] != "between"
                                  # 2026-08-12（系统性审查，P2）：float() 幂等防御——
                                  # _load_rules 加载时已归一化，此处与 _eval_rule 同口径
                                  else f"[{float(r['low'])}, {float(r['high'])})")
        else:
            # 2026-08-13（二审 #8）：单阈值趋势补 ">= " 前缀——与 _eval_rule 的
            # threshold 展示（">= 50"）一致；此前输出 "up 50%" 缺比较符，语义模糊。
            entry["criterion"] = (f"{r['direction']} >= {r['threshold_pct']}%"
                                  if "low_pct" not in r
                                  else f"{r['direction']} [{r['low_pct']}, {r['high_pct']})%")
        out.append(entry)
    return {"ok": True, "data": {"rules": out}}


def explain_verdict(evaluation: Dict[str, Any]) -> Dict[str, Any]:
    """把 evaluate_risk_rules 的结果翻成判定链路，供审计与医生复核。

    入参接受 evaluate_risk_rules 的 {ok,data} 信封或直接其 data 体（向后兼容）。
    BUG-58（2026-08-12）：显式失败信封（{ok:false}）与非 dict 入参直接拒绝——
    此前失败信封会被误读为"无规则命中"，掩盖评估失败的事实。
    """
    # 2026-08-12（系统性审查）：补鉴权——本函数与 list_rules 是模块内仅有的两个
    # 未走 enforce_read 的对外函数（server 工具层同样无 enforce），矩阵收紧将失效。
    caller = get_caller()
    enforce_read(MCP_NAME)
    if not isinstance(evaluation, dict):
        raise ValueError(f"无法解释无效的评估结果：{evaluation!r}（需要评估结果 dict）")
    if evaluation.get("ok") is False:
        raise ValueError(
            f"无法解释失败的评估结果：{evaluation.get('error')} {evaluation.get('detail', '')}".strip())
    # 信封模式（含 data 键）：data 必须为 dict——{ok:true, data:null} 属于异常信封，
    # 显式报错而非静默当"无规则命中"（fail-closed，与 BUG-58 原则一致）。
    # 裸 data 模式（无 data 键）：evaluation 即业务体（向后兼容）。
    if "data" in evaluation:
        data = evaluation["data"]
        if not isinstance(data, dict):
            raise ValueError(
                f"无法解释无效的评估结果：{evaluation.get('ok')} 信封的 data 应为 dict，"
                f"实际为 {type(data).__name__}")
        payload = data
    else:
        payload = evaluation
    chain = []
    # BUG（2026-08-12，P0）：assess_clinical_status 输出的命中规则键名为
    # risk_matched_rules，而 evaluate_risk_rules 用 matched_rules——直接把 DAG 结果
    # 传给 explain_verdict 会静默回退"无规则命中"，高危病人被误判"未触发任何规则"。
    # 双键名兼容（元素结构相同：均来自 evaluate_risk_rules 的 matched_rules）。
    # 2026-08-12（系统性审查，P3）：显式按键存在性提取——此前 `get("matched_rules") or
    # get("risk_matched_rules", [])` 依赖 Python Falsy：matched_rules 为 []（已评估但
    # 无命中）时也会强算右侧表达式，属隐式类型评估；现保留"键存在且为 [] = 明确无命中"
    # 语义，仅当键缺失或为 None 时回退另一键名（[] 不再被 or 吞掉）。
    matched_rules = payload.get("matched_rules")
    if matched_rules is None:
        matched_rules = payload.get("risk_matched_rules", [])
    # 2026-08-12（系统性审查，P2）：非法 evaluation 的结构校验——LLM/客户端可能传入
    # 残缺/错误类型的 matched_rules（如 [{"id":"R-01"}] 缺 name、元素为 str、整体非 list）。
    # 此前直接索引 m["id"]/m["name"] 抛 KeyError/TypeError → _invalid 归 INTERNAL_ERROR，
    # 客户端输入错误触发服务端误报监控。现显式结构校验：非法 → ValueError（→ INVALID_INPUT），
    # 与"内部代码 bug（键名写错）仍抛 KeyError → INTERNAL_ERROR"的监控语义保持区分。
    if not isinstance(matched_rules, list):
        raise ValueError(
            f"无法解释无效的评估结果：matched_rules 应为列表，实际为 {type(matched_rules).__name__}"
        )
    _MATCHED_RULE_FIELDS = ("id", "name", "level", "observed", "threshold", "unit", "description")
    for i, m in enumerate(matched_rules):
        if not isinstance(m, dict):
            raise ValueError(f"matched_rules[{i}] 应为 dict，实际为 {type(m).__name__}")
        missing = [f for f in _MATCHED_RULE_FIELDS if f not in m]
        if missing:
            raise ValueError(f"matched_rules[{i}] 缺少字段: {missing}")
    for m in matched_rules:
        chain.append({
            "rule_id": m["id"],
            "rule_name": m["name"],
            "level": m["level"],
            "observed": m["observed"],
            "threshold": m["threshold"],
            "unit": m["unit"],
            "why": m["description"],
        })
    if not chain:
        chain.append({"rule_id": "-", "rule_name": "无规则命中",
                      "level": "none", "observed": None, "threshold": None,
                      "unit": None, "why": "本轮数据未触发任何风险规则"})
    return {"ok": True, "data": {"chain": chain}}

# ================================================================
# DAG: assess_clinical_status (v2.2/v2.3)
# ================================================================

def assess_clinical_status(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    serum_creatinine_unit: str = "mg_dL",
    is_preterm: bool = False,
    sex: Optional[str] = None,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    method: Optional[EgfrMethod] = None,
    new_labs: Optional[Dict[str, float]] = None,
    prior_labs: Optional[Dict[str, float]] = None,
    prior_level: Optional[str] = None,
    dialysis_mode: Optional[str] = None,
    upcr_mg_mmol: Optional[float] = None,
) -> Dict[str, Any]:
    """一键评估 CKD 临床状态（eGFR + 分期 + 风险评分 DAG）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1）。
    method=None 时自动推理（M-1 后：有 bun → 床旁式+告警，有 k_value → classic，
    <1y 早产儿 → classic[k=0.33]，否则 bedside2009）；传入 method 则优先使用传入值。
    sex（M/F/male/female/男/女）仅 method="classic" 且 ≥13y 时生效（女性 k=0.55）。
    serum_creatinine_unit（BUG-40 修复）：mg_dL 默认 / umol_L 自动 ÷88.4——P1 get_labs
    返回 scr_umol_L（µmol/L），直接透传会导致 eGFR 与风险规则的 scr 判断同时错 88 倍。
    dialysis_mode（S1 修复，2026-08-13）：none/hemodialysis/peritoneal——透析患儿
    eGFR<15 分期为 G5D，供 care 层随访频率区分（G5D 每月 vs G5 每 60 天）。
    BUG-34：显式取 caller 回写审计字段。
    BUG-35 说明（2026-08-12）：DAG 入口 enforce_read 后，内部子函数（calc_egfr_schwartz /
    classify_ckd / evaluate_risk_rules）各自再 enforce_read——这是**有意的防御纵深**
    （子函数可被外部直接调用，不能假设必经 DAG），重复校验为 O(1) 查矩阵，开销可忽略。
    """
    caller = get_caller()  # BUG-34：显式取 caller 回写审计字段
    enforce_read(MCP_NAME)
    # F5（2026-08-16，第七轮审查）：DAG 入口 dialysis_mode 白名单（防御纵深）——
    # 正常路径经 classify_ckd → _egfr_to_g 叶子校验，但 DAG 重构绕过叶子即漏；
    # 入口同口径显式拒绝非法值（与 _egfr_to_g 一致）。
    if dialysis_mode is not None:
        if not isinstance(dialysis_mode, str) \
                or dialysis_mode.strip().lower() not in ("none", "hemodialysis", "peritoneal"):
            raise ValueError(
                f"dialysis_mode 必须是 none / hemodialysis / peritoneal 之一，收到："
                f"{dialysis_mode!r}")

    # 2026-08-12（系统性审查，P1）：sex 入口统一清洗为 "male"/"female"/None——
    # 与 calc 内部归一化幂等，保证全链路口径一致（审计/日志/sex_note 单一事实源）。
    sex = _normalize_sex(sex)

    # 2026-08-12（系统性审查，P2）：弱必填形参 Fail-Fast 强校验——uacr/upcr/bun/k_value
    # 此前在 labs 构建中途才被 _require_finite 拦截，非法值（NaN/Inf）会带着进入中间
    # 计算（method 自动推理等）才报错；现于入口统一校验（与 age/height 同口径），
    # 非法输入即刻失败（INVALID_INPUT），不进入任何计算。scr/egfr 由 calc 内部保证有限。
    if uacr_mg_g is not None:
        uacr_mg_g = _require_finite(uacr_mg_g, "uacr_mg_g")
    if upcr_mg_g is not None:
        upcr_mg_g = _require_finite(upcr_mg_g, "upcr_mg_g")
    if bun_mg_dl is not None:
        bun_mg_dl = _require_finite(bun_mg_dl, "bun_mg_dl")
    if k_value is not None:
        k_value = _require_finite(k_value, "k_value")

    # BUG-40：DAG 内统一归一化肌酐到 mg/dL（eGFR 公式与风险规则 scr 判断共用同一值）
    scr_mgdl = _normalize_scr(serum_creatinine_mgdl, serum_creatinine_unit)

    # 2026-08-12（系统性审查，P2）：显式形参与 new_labs 字典的双向对齐——调用方若只在
    # new_labs 传 {"bun"/"uacr"/"upcr"} 而未传形参，此前会导致：① method 自动推理看不到
    # BUN → 误降级 bedside2009；② classify_ckd 收不到 uacr/upcr → 白蛋白尿 A 分期丢失。
    # 形参为 None 时从归一化 new_labs 补齐（_normalize_labs 已做有限性强校验，脏值在此抛错）。
    norm_inputs = _normalize_labs(new_labs) or {}
    if bun_mg_dl is None and "bun" in norm_inputs:
        bun_mg_dl = norm_inputs["bun"]
    if uacr_mg_g is None and "uacr" in norm_inputs:
        uacr_mg_g = norm_inputs["uacr"]
    if upcr_mg_g is None and "upcr" in norm_inputs:
        upcr_mg_g = norm_inputs["upcr"]

    # 自动识别 Schwartz 方法（显式传入优先）
    warnings: List[str] = []
    if method is None:
        # M-1（2026-08-15）：提供 BUN 不再自动默认 revised2009（假公式已移除）——
        # 含 BUN 修订需胱抑素 C（CKiD 组合式）本系统未实现，BUN 不参与计算，降级床旁式
        # 并显式告警（不静默吞参）。
        if bun_mg_dl is not None:
            method = "bedside2009"
            warnings.append(
                "提供了 bun_mg_dl 但未指定 method：含 BUN 的修订需胱抑素 C（CKiD "
                "组合式），本系统未实现，BUN 不参与计算；已按床旁式（0.413）处理。"
                "如需自定义 k 值请显式 method='classic' 并传 k_value。")
            if k_value is not None:
                warnings.append(
                    "同时提供 bun_mg_dl 与 k_value 且未指定 method：k_value 被忽略；"
                    "如需自定义 k 值请显式 method='classic' 并传 k_value。")
        elif k_value is not None:
            method = "classic"
        elif is_preterm and age_years < 1:
            # 2026-08-12（系统性审查，P1）：早产儿自动推理修正——此前 is_preterm=True
            # 且无 bun/k_value 时落到 bedside2009（固定 k=0.413），早产 k=0.33 静默失效，
            # eGFR 高估约 25%（0.413/0.33），G4 可能误判 G3，违背 BUG-47 防护初衷。
            # <1y 早产儿自动切经典式（k=0.33 生效），与 BUG-47 口径一致。
            method = "classic"
            warnings.append("检测到早产（is_preterm=True 且 <1 岁）：自动采用经典式 "
                            "k=0.33；如需床旁式请显式 method='bedside2009'。")
        else:
            method = "bedside2009"

    egfr_r = calc_egfr_schwartz(
        age_years=age_years,
        height_cm=height_cm,
        # 2026-08-12（系统性审查，P2）：透传**原始**肌酐值与单位——此前传预归一化的
        # scr_mgdl + 默认 unit="mg_dL"，calc 返回的审计字段 scr_unit_used 恒为 "mg/dL"，
        # 用户原始输入 umol_L 的追溯信息被抹除。calc 内部会自行归一化（_normalize_scr），
        # labs["scr"] 仍用下方归一化值，两者互不影响。
        serum_creatinine_mgdl=serum_creatinine_mgdl,
        serum_creatinine_unit=serum_creatinine_unit,
        method=method,
        bun_mg_dl=bun_mg_dl,
        k_value=k_value,
        is_preterm=is_preterm,
        sex=sex,
    )
    # BUG-67 后补（2026-08-12）：防御性兜底——calc_egfr_schwartz 当前对非法输入抛 ValueError
    # （冒泡到 server._invalid 归 INVALID_INPUT），但 DAG 内部直接索引 egfr_r["data"] 依赖
    # "成功必返回 dict 信封"的隐含假设；若未来底层重构为"错误返回错误信封"而非抛异常，
    # 此处会 KeyError。统一先判 ok 再透传，双形态兼容。
    if not egfr_r.get("ok"):
        return egfr_r
    ckd_r = classify_ckd(
        # N1 修复（2026-08-13）：判级用未舍入原始值 egfr_raw——round 值在边界
        # （如 89.96→90.0）会翻 G2→G1；展示仍用 egfr（round 0.1）。
        egfr=egfr_r["data"].get("egfr_raw", egfr_r["data"]["egfr"]),
        uacr_mg_g=uacr_mg_g,
        upcr_mg_g=upcr_mg_g,
        dialysis_mode=dialysis_mode,
        upcr_mg_mmol=upcr_mg_mmol,  # M-5（2026-08-16）：UPCR mg/mmol 入参透传
    )
    # 始终以本轮计算的 eGFR + 已传入参数打底做风险评估，不因未传 new_labs 而静默跳过
    # 2026-08-12（系统性审查，P3）：labs 直接基于**已归一化**的 norm_inputs 构造——
    # 此前基于原始 new_labs 再调 _normalize_labs(labs) 重复一遍完整键名换算+遍历。
    # norm_inputs 已含完整键名→短名转换与 _require_finite 强校验；此处仅并入 DAG
    # 计算的短名键（scr/egfr，calc 已保证有限）与形参键（bun/uacr/upcr，入口已
    # Fail-Fast 校验为有限值或 None；形参优先于 new_labs 值），无需重复校验。
    labs = dict(norm_inputs)
    # N1 修复：R-07 趋势规则同样用未舍入 egfr_raw 判级（round 值会扰动 -25% 阈值判定）
    _egfr_for_rules = egfr_r["data"].get("egfr_raw", egfr_r["data"]["egfr"])
    for k, v in (("scr", scr_mgdl), ("egfr", _egfr_for_rules),
                 ("bun", bun_mg_dl), ("uacr", uacr_mg_g), ("upcr", upcr_mg_g)):
        if v is not None:
            labs[k] = v
    # 2026-08-12（系统性审查，P3）：prior_labs 归一化前移一次，evaluate 与趋势完整性
    # 校验共用（evaluate 内部再归一化幂等无副作用）——此前两处各归一化一遍。
    norm_prior = _normalize_labs(prior_labs)
    risk_r = evaluate_risk_rules(
        new_labs=labs,
        prior_labs=norm_prior,
        prior_level=prior_level,
    )
    risk_data = risk_r.get("data", {}) if risk_r.get("ok") else {}
    # BUG-26 修复（2026-08-12）：显式声明本轮风险扫描的覆盖完整性——
    # rules.json 全部规则依赖 ca/egfr/hb/k/p/scr/ua，若调用方未传 new_labs，
    # labs 仅含 scr+egfr（+可选 bun/uacr/upcr），电解质/贫血规则会被静默跳过，
    # 直接返回 overall_level="none" 会给临床"已全面评估且无风险"的假象。
    # 规则库读一次复用（避免 DAG 内多路各 deepcopy 一次）。
    rules_doc = _load_rules()
    rule_metrics = sorted({rule["metric"] for rule in rules_doc["rules"]})
    # BUG-58（2026-08-12）：覆盖度统计须基于归一化键名——调用方传 P1 完整键名
    # （k_mmol_L、hb_g_L 等）时，原始 labs 不含短名，直接求差集会误报"未覆盖
    # K/Hb"，与实际已命中的规则矛盾（用户审查发现）。
    # 2026-08-12（系统性审查，P3）：labs 已基于 norm_inputs 构造（归一化完成），
    # 直接复用，不再重复 _normalize_labs(labs)（省一次完整遍历）。
    normalized_labs = labs
    missing_metrics = sorted(set(rule_metrics) - set(normalized_labs))
    covered_metrics = sorted(set(rule_metrics) - set(missing_metrics))
    if missing_metrics:
        risk_completeness = {
            "covered_metrics": covered_metrics,
            "missing_metrics": missing_metrics,
            "note": (f"本轮仅评估了 {covered_metrics} 相关规则，"
                     f"未覆盖指标: {missing_metrics}；如涉及电解质（K/Ca/P）或贫血（Hb）危急值，"
                     "请传入 new_labs 补充后重评，overall_level=none 不代表全面无风险。"),
        }
    else:
        risk_completeness = {
            "covered_metrics": rule_metrics,
            "missing_metrics": [],
            "note": "本轮已覆盖规则引擎全部依赖指标。",
        }
    # 2026-08-12（系统性审查）：覆盖度假象修正——rules.json 含 trend_pct 规则
    # （R-01/R-08 scr 环比、R-07 egfr 环比），它们依赖 prior_labs 历史对照；
    # 若未传 prior_labs，missing_metrics 即使为空，"已全面覆盖"也是假象
    # （趋势规则被 _eval_rule 静默跳过）。显式提示，不给临床"已全面评估"的承诺。
    # 2026-08-12（系统性审查，P2 修正上轮粗判）：趋势规则完整性须**精确校验其依赖的
    # 具体指标**在 prior_labs 中是否覆盖——此前 len(prior_labs)>0 只判"非空"：调用方传
    # prior_labs={"scr": 1.0} 时 has_prior=True，但 R-07 依赖 egfr（prior 缺 egfr）仍被
    # 静默跳过且无提示。现按趋势规则 metric 集（scr/egfr）逐一核对归一化 prior_labs。
    trend_rules = [r for r in rules_doc["rules"] if r["type"] == "trend_pct"]
    if trend_rules:
        required_prior = {r["metric"] for r in trend_rules}
        norm_prior = norm_prior or {}  # 复用前移的归一化结果（评估与完整性校验共用）
        missing_prior = sorted(required_prior - set(norm_prior))
        if missing_prior:
            risk_completeness["note"] += (
                f"（注意：prior_labs 缺失 {missing_prior} 历史对照，对应动态趋势类规则"
                "（R-01/R-07/R-08）未触发评估，overall_level 不含相关急性恶化判断。）")
    return {
        "ok": True,
        "data": {
            "caller": caller,  # BUG-34 审计字段
            "egfr": egfr_r["data"]["egfr"],
            "egfr_unit": egfr_r["data"]["unit"],
            "egfr_method": egfr_r["data"]["method"],
            # 四审（2026-08-12）：透出计算层警示——calc 的 note（早产 k 覆盖提示、
            # BUN 静默丢弃提示、经典式参数说明等）与 classify 的 risk_note（进展风险
            # 提示）此前仅在底层函数返回，DAG 结果里不可见；医生看一键评估结论会漏掉
            # 关键口径警示（如"is_preterm 用床旁式无早产修正"）。
            "egfr_note": egfr_r["data"].get("note"),
            "ckd_risk_note": ckd_r["data"].get("risk_note"),
            "ckd_stage": ckd_r["data"]["stage"],
            "ckd_g_stage": ckd_r["data"]["g"],
            "ckd_a_stage": ckd_r["data"].get("a"),
            "risk_level": risk_data.get("overall_level", "none"),
            "risk_matched_rules": risk_data.get("matched_rules", []),
            # BUG-21 修复：DAG 透出历史等级对比（prior_comparison），此前被丢弃
            "prior_comparison": risk_data.get("prior_comparison"),
            # BUG-26 修复：风险扫描覆盖完整性声明（防"未全面评估却显示无风险"的假象）
            "risk_completeness": risk_completeness,
            # 2026-08-12：DAG 级警告（如 k_value 被忽略）
            "warnings": warnings,
        },
    }

# ---- M8: risk rules engine ----
