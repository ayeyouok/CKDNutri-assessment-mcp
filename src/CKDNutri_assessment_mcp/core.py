# -*- coding: utf-8 -*-
"""M6 临床计算纯函数：eGFR-Schwartz 与 KDIGO 儿童 CKD 分期。

不依赖 fastmcp，可直接 import 单测。所有公式口径在 ADR-006 锁定。

公式来源与口径（ADR-006）：
- eGFR（床旁 Schwartz 2009）：eGFR = 0.413 × 身高(cm) / Scr(mg/dL)，结果以 ml/min/1.73m² 计。
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
import os

from typing import Any, Dict, List, Literal, Optional

from a207_policy import enforce_read, get_caller

MCP_NAME = "CKDNutri-assessment-mcp"

EgfrMethod = Literal["bedside2009", "revised2009", "classic"]

_BEDSIDE_K = 0.413  # Schwartz 2009 床旁系数


def _egfr_to_g(egfr: float) -> str:
    """KDIGO 2024 儿童 eGFR 分期（ml/min/1.73m²）。"""
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


def _schwartz_k(age_years: float, k_value: Optional[float], is_preterm: bool = False) -> float:
    """Schwartz 经典 k 值：早产儿(<37周) 0.33、<1yr 足月儿 0.45、1-12yr=0.55、≥13yr=0.70。

    BUG-47（2026-08-12）：早产儿 k=0.33 正式实现——CAKUT（先天性肾发育不良）是儿童 CKD
    首要病因，早产儿占比高；用 0.45 替代 0.33 会把 eGFR 高估约 36%（0.45/0.33），
    可能将 G4 误判为 G3。调用方须在已知早产时传 is_preterm=True（<1yr 生效）。
    """
    if k_value is not None:
        return float(k_value)
    if age_years < 1:
        return 0.33 if is_preterm else 0.45
    if age_years < 13:
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
) -> dict:
    """估算肾小球滤过率（Schwartz 系列）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    返回：egfr(ml/min/1.73m²)、method、formula(所用公式串)、note(警示/口径)。
    BUG-34（2026-08-12）：显式取 caller 并回写审计字段（此前仅 enforce_read 内部取用，
    调用者身份不落返回，无法追溯谁触发了计算）。
    BUG-40（2026-08-12）：新增 serum_creatinine_unit 单位归一化——P1 clinical-data
    `get_labs` 返回 `scr_umol_L`（µmol/L），本函数历史上只接受 mg/dL；编排层若把 P1
    数值直接透传会错 88.4 倍（eGFR 缩小 88 倍）。现支持 unit="umol_L" 自动 ÷88.4 转
    mg/dL（默认仍为 mg_dL，向后兼容）。也提供显式转换函数 `scr_umol_to_mgdl`。
    BUG-47（2026-08-12）：新增 is_preterm——早产儿经典 k=0.33（仅 <1yr 生效），
    防 eGFR 高估 36%。
    """
    caller = get_caller()
    enforce_read(MCP_NAME)
    if method not in ("bedside2009", "revised2009", "classic"):
        raise ValueError(
            f"无效的 method: '{method}'，必须为 bedside2009 / revised2009 / classic 之一"
        )
    if age_years < 0:
        raise ValueError("age_years 不能为负")
    if height_cm <= 0:
        raise ValueError("height_cm 必须 > 0")
    scr = _normalize_scr(serum_creatinine_mgdl, serum_creatinine_unit)
    if scr <= 0:
        raise ValueError("serum_creatinine 必须 > 0")
    if method == "revised2009" and (bun_mg_dl is None or bun_mg_dl <= 0):
        raise ValueError("revised2009 需要 bun_mg_dl > 0")

    if method == "classic":
        k = _schwartz_k(age_years, k_value, is_preterm=is_preterm)
        egfr = k * height_cm / scr
        formula = f"eGFR = k×height/Scr, k={k}"
        k_note = "早产儿 k=0.33" if (is_preterm and age_years < 1) else "经典年龄带 k"
        note = f"{k_note}；<1y=0.45, 1-12y=0.55, ≥13y=0.70，可被 k_value 覆盖。"
    elif method == "revised2009":
        # bun_mg_dl 已由前置校验保证非空 > 0，此处不再重复判定
        denom = scr + 0.003 * bun_mg_dl - 0.024
        if denom <= 0:
            raise ValueError("revised2009 分母非正（Scr+0.003×BUN-0.024 必须 > 0）")
        egfr = _BEDSIDE_K * height_cm / denom
        formula = f"eGFR = 0.413×height/(Scr + 0.003×BUN − 0.024)"
        note = "含 BUN 修订 Schwartz 2009，对 eGFR<60 的儿童更准确。"
    else:  # bedside2009
        egfr = _BEDSIDE_K * height_cm / scr
        formula = "eGFR = 0.413×height/Scr"
        note = "床旁 Schwartz 2009（KDIGO 推荐默认式）。"

    egfr = round(egfr, 1)
    pediatric_caveat = (
        "注意：<2 岁婴儿 eGFR 参考范围低于年长儿（正常可低至 60–90），"
        "G1/G2 阈值在婴儿期需结合月龄与生长曲线解读。"
        if age_years < 2 else ""
    )
    # BUG-15：成功返回统一 {ok, data} 信封（此前扁平结构）；BUG-34：含 caller 审计字段
    return {
        "ok": True,
        "data": {
            "caller": caller,
            "egfr": egfr,
            "unit": "ml/min/1.73m2",
            "method": method,
            "formula": formula,
            "pediatric_caveat": pediatric_caveat,
            "note": note,
            "scr_unit_used": _unit_label(serum_creatinine_unit),
        },
    }


def scr_umol_to_mgdl(value_umol_L: float) -> float:
    """显式单位转换：µmol/L → mg/dL（÷88.4）。供编排层把 P1 get_labs 的
    `scr_umol_L` 转成 P4 eGFR 公式所需单位（BUG-40）。"""
    return float(value_umol_L) / 88.4


def _normalize_scr(value: float, unit: str) -> float:
    """按声明单位归一化肌酐到 mg/dL；非法单位显式报错（fail-closed）。"""
    v = float(value)
    u = (unit or "mg_dL").strip().lower().replace("µ", "u")
    if u in ("mg_dl", "mgdl", "mg"):
        return v
    if u in ("umol_l", "umoll", "umol"):
        return v / 88.4
    raise ValueError(f"无效的 serum_creatinine_unit: {unit!r}，可用 mg_dL / umol_L")


def _unit_label(unit: str) -> str:
    u = (unit or "mg_dL").strip().lower().replace("µ", "u")
    return "umol/L（已 ÷88.4 转 mg/dL）" if u in ("umol_l", "umoll", "umol") else "mg/dL"


def classify_ckd(
    egfr: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
) -> dict:
    """KDIGO 2024 儿童 CKD 合并分期。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    egfr 必需；白蛋白尿二选一（uacr 或 upcr）。返回 g、a、stage(GxAx)、description、risk_note。
    BUG-34：显式取 caller 并回写审计字段。
    """
    caller = get_caller()
    enforce_read(MCP_NAME)
    if egfr < 0:
        raise ValueError("egfr 必须 >= 0")
    if uacr_mg_g is not None and uacr_mg_g < 0:
        raise ValueError("uacr_mg_g 不能为负")
    if upcr_mg_g is not None and upcr_mg_g < 0:
        raise ValueError("upcr_mg_g 不能为负")

    g = _egfr_to_g(egfr)
    a: Optional[str] = None
    albuminuria_source = ""
    albuminuria_note = None
    if uacr_mg_g is not None:
        a = _acr_to_a(uacr_mg_g)
        albuminuria_source = "UACR"
    elif upcr_mg_g is not None:
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
    }[g]
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
            "a_description": a_desc.get(a) if a else None,
            "albuminuria_source": albuminuria_source or None,
            "albuminuria_note": albuminuria_note,
            "risk_note": risk_note,
        },
    }


def _risk_note(g: str, a: Optional[str]) -> str:
    """按 G/A 给出进展风险与随访强度提示（信息性）。"""
    g_rank = {"G1": 0, "G2": 1, "G3a": 2, "G3b": 3, "G4": 4, "G5": 5}[g]
    a_rank = {"A1": 0, "A2": 1, "A3": 2}.get(a or "A1", 0)
    score = g_rank + a_rank
    if score >= 6:
        return "进展风险高：建议缩短随访间隔（如 1–3 个月）并由肾科密切管理。"
    if score >= 3:
        return "进展风险中等：建议 3–6 个月随访一次。"
    return "进展风险低：建议 6–12 个月常规随访。"
_LEVEL_RANK = {"L1": 3, "L2": 2, "L3": 1, "none": 0}

_RULES_PATH = os.path.join(os.path.dirname(__file__), "data", "rules.json")
_RULES: Optional[Dict[str, Any]] = None


def _load_rules() -> Dict[str, Any]:
    global _RULES
    if _RULES is None:
        try:
            with open(_RULES_PATH, "r", encoding="utf-8") as f:
                _RULES = json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"风险规则文件缺失：{_RULES_PATH}；请确认 data/rules.json 存在"
            )
        except json.JSONDecodeError as e:
            raise ValueError(
                f"风险规则文件 JSON 解析失败：{_RULES_PATH}，{e}"
            )
    return _RULES


def _pct_change(new: float, old: float) -> float:
    """相对变化百分比（正=升，负=降）。old<=0 视为无效。"""
    try:
        new, old = float(new), float(old)
    except (ValueError, TypeError):
        return float("nan")
    if old is None or old <= 0:
        return float("nan")
    return (new - old) / old * 100.0


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
            hit = False
        if not hit:
            return None
        op_label = {"gt": ">", "gte": ">=", "lt": "<", "lte": "<="}.get(op, op)
        return {
            "metric": metric,
            "observed": v,
            "threshold": (rule.get("low"), rule.get("high")) if op == "between"
                         else f"{op_label} {rule['threshold']}",
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
            "observed": round(pct, 1),
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
    "scr_umol_L": ("scr", 1.0 / 88.4),   # µmol/L → mg/dL
    "k_mmol_L": ("k", 1.0),
    "p_mmol_L": ("p", 1.0),
    "hb_g_L": ("hb", 1.0),
    "ca_mmol_L": ("ca", 1.0),
    "na_mmol_L": ("na", 1.0),
    "ua_umol_L": ("ua", 1.0),
    "egfr_ml_min": ("egfr", 1.0),
    "bun_mmol_L": ("bun", 1.0),
}


def _normalize_labs(labs: Optional[Dict[str, float]]) -> Optional[Dict[str, float]]:
    """把 P1 完整键名 + PCP 单位归一化为规则短名 + 规则单位（BUG-45）。

    已用短名的输入原样保留；完整键名在短名缺失时转换并换算单位。
    返回新 dict，不修改调用方对象。
    """
    if not labs:
        return labs
    out: Dict[str, float] = dict(labs)
    for full_key, (short, factor) in _LAB_ALIAS_TO_RULE.items():
        if full_key in out and short not in out:
            try:
                out[short] = float(out[full_key]) * factor
            except (TypeError, ValueError):
                continue
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

    # 对比：仅展示，不兜底
    delta_note = "无历史等级对照"
    if prior_level:
        if prior_level == overall:
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
            "prior_comparison": {
                "prior_level": prior_level,
                "current_level": overall,
                "delta_note": delta_note,
            },
            "level_correction_applied": True,
        },
    }


def list_rules() -> List[Dict[str, Any]]:
    """返回规则清单（不含评估逻辑）。"""
    rules_doc = _load_rules()
    out = []
    for r in rules_doc["rules"]:
        entry = {"id": r["id"], "name": r["name"], "level": r["level"],
                 "type": r["type"], "description": r["description"], "unit": r["unit"]}
        if r["type"] == "absolute":
            entry["criterion"] = (f"{r['operator']} {r.get('threshold')}"
                                  if r["operator"] != "between"
                                  else f"[{r['low']}, {r['high']})")
        else:
            entry["criterion"] = (f"{r['direction']} {r.get('threshold_pct', '')}%"
                                  if "low_pct" not in r
                                  else f"{r['direction']} [{r['low_pct']}, {r['high_pct']})%")
        out.append(entry)
    return {"ok": True, "data": {"rules": out}}


def explain_verdict(evaluation: Dict[str, Any]) -> Dict[str, Any]:
    """把 evaluate_risk_rules 的结果翻成判定链路，供审计与医生复核。

    入参接受 evaluate_risk_rules 的 {ok,data} 信封或直接其 data 体（向后兼容）。
    """
    data = evaluation.get("data") if isinstance(evaluation, dict) else None
    payload = data if isinstance(data, dict) else evaluation
    chain = []
    for m in payload.get("matched_rules", []):
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
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    method: Optional[EgfrMethod] = None,
    new_labs: Optional[Dict[str, float]] = None,
    prior_labs: Optional[Dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> Dict[str, Any]:
    """一键评估 CKD 临床状态（eGFR + 分期 + 风险评分 DAG）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1）。
    method=None 时自动推理（有 bun → revised2009，有 k_value → classic，否则 bedside2009）；
    传入 method 则优先使用传入值。
    serum_creatinine_unit（BUG-40 修复）：mg_dL 默认 / umol_L 自动 ÷88.4——P1 get_labs
    返回 scr_umol_L（µmol/L），直接透传会导致 eGFR 与风险规则的 scr 判断同时错 88 倍。
    BUG-34：显式取 caller 回写审计字段。
    BUG-35 说明（2026-08-12）：DAG 入口 enforce_read 后，内部子函数（calc_egfr_schwartz /
    classify_ckd / evaluate_risk_rules）各自再 enforce_read——这是**有意的防御纵深**
    （子函数可被外部直接调用，不能假设必经 DAG），重复校验为 O(1) 查矩阵，开销可忽略。
    """
    caller = get_caller()  # BUG-34：显式取 caller 回写审计字段
    enforce_read(MCP_NAME)

    # BUG-40：DAG 内统一归一化肌酐到 mg/dL（eGFR 公式与风险规则 scr 判断共用同一值）
    scr_mgdl = _normalize_scr(serum_creatinine_mgdl, serum_creatinine_unit)

    # 自动识别 Schwartz 方法（显式传入优先）
    if method is None:
        if bun_mg_dl is not None:
            method = "revised2009"
        elif k_value is not None:
            method = "classic"
        else:
            method = "bedside2009"

    egfr_r = calc_egfr_schwartz(
        age_years=age_years,
        height_cm=height_cm,
        serum_creatinine_mgdl=scr_mgdl,
        method=method,
        bun_mg_dl=bun_mg_dl,
        k_value=k_value,
        is_preterm=is_preterm,
    )
    ckd_r = classify_ckd(
        egfr=egfr_r["data"]["egfr"],
        uacr_mg_g=uacr_mg_g,
        upcr_mg_g=upcr_mg_g,
    )
    risk_r: Dict[str, Any] = {}
    # 始终以本轮计算的 eGFR + 已传入参数打底做风险评估，不因未传 new_labs 而静默跳过
    labs = dict(new_labs) if new_labs else {}
    labs["scr"] = scr_mgdl
    labs["egfr"] = egfr_r["data"]["egfr"]
    if bun_mg_dl is not None:
        labs["bun"] = bun_mg_dl
    if uacr_mg_g is not None:
        labs["uacr"] = uacr_mg_g
    if upcr_mg_g is not None:
        labs["upcr"] = upcr_mg_g
    risk_r = evaluate_risk_rules(
        new_labs=labs,
        prior_labs=prior_labs,
        prior_level=prior_level,
    )
    risk_data = risk_r.get("data", {}) if risk_r.get("ok") else {}
    # BUG-26 修复（2026-08-12）：显式声明本轮风险扫描的覆盖完整性——
    # rules.json 全部规则依赖 ca/egfr/hb/k/p/scr/ua，若调用方未传 new_labs，
    # labs 仅含 scr+egfr（+可选 bun/uacr/upcr），电解质/贫血规则会被静默跳过，
    # 直接返回 overall_level="none" 会给临床"已全面评估且无风险"的假象。
    rule_metrics = sorted({rule["metric"] for rule in _load_rules()["rules"]})
    missing_metrics = sorted(set(rule_metrics) - set(labs))
    if missing_metrics:
        risk_completeness = {
            "covered_metrics": sorted(set(rule_metrics) - set(missing_metrics)),
            "missing_metrics": missing_metrics,
            "note": (f"本轮仅评估了 {sorted(set(rule_metrics) - set(missing_metrics))} 相关规则，"
                     f"未覆盖指标: {missing_metrics}；如涉及电解质（K/Ca/P）或贫血（Hb）危急值，"
                     "请传入 new_labs 补充后重评，overall_level=none 不代表全面无风险。"),
        }
    else:
        risk_completeness = {
            "covered_metrics": rule_metrics,
            "missing_metrics": [],
            "note": "本轮已覆盖规则引擎全部依赖指标。",
        }
    return {
        "ok": True,
        "data": {
            "caller": caller,  # BUG-34 审计字段
            "egfr": egfr_r["data"]["egfr"],
            "egfr_unit": egfr_r["data"]["unit"],
            "egfr_method": egfr_r["data"]["method"],
            "ckd_stage": ckd_r["data"]["stage"],
            "ckd_g_stage": ckd_r["data"]["g"],
            "ckd_a_stage": ckd_r["data"].get("a"),
            "risk_level": risk_data.get("overall_level", "none"),
            "risk_matched_rules": risk_data.get("matched_rules", []),
            # BUG-21 修复：DAG 透出历史等级对比（prior_comparison），此前被丢弃
            "prior_comparison": risk_data.get("prior_comparison"),
            # BUG-26 修复：风险扫描覆盖完整性声明（防"未全面评估却显示无风险"的假象）
            "risk_completeness": risk_completeness,
        },
    }

# ---- M8: risk rules engine ----
