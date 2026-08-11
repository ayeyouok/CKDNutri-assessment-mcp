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


def _pcr_to_a(upcr_mg_g: float) -> str:
    """蛋白尿 A 期：UPCR mg/g（儿童常以 PCR 评估）。"""
    if upcr_mg_g < 150:
        return "A1"
    if upcr_mg_g <= 500:
        return "A2"
    return "A3"


def _schwartz_k(age_years: float, k_value: Optional[float]) -> float:
    """Schwartz 经典 k 值：<1yr=0.45, 1-12yr=0.55, ≥13yr=0.70。
    注：早产儿（<37周）理论上需 0.33，本函数暂不内置早产判定（需调用方传入 is_preterm
    并自行调整 k_value），<1yr 默认使用足月儿 0.45。
    """
    if k_value is not None:
        return float(k_value)
    if age_years < 1:
        return 0.45
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
) -> dict:
    """估算肾小球滤过率（Schwartz 系列）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    返回：egfr(ml/min/1.73m²)、method、formula(所用公式串)、note(警示/口径)。
    """
    enforce_read(MCP_NAME)
    if age_years < 0:
        raise ValueError("age_years 不能为负")
    if height_cm <= 0:
        raise ValueError("height_cm 必须 > 0")
    if serum_creatinine_mgdl <= 0:
        raise ValueError("serum_creatinine_mgdl 必须 > 0")
    if method == "revised2009" and (bun_mg_dl is None or bun_mg_dl <= 0):
        raise ValueError("revised2009 需要 bun_mg_dl > 0")

    if method == "classic":
        k = _schwartz_k(age_years, k_value)
        egfr = k * height_cm / serum_creatinine_mgdl
        formula = f"eGFR = k×height/Scr, k={k}"
        note = "经典 k 值 Schwartz；k 默认按年龄带（<1y=0.45, 1-12y=0.55, ≥13y=0.70），可被 k_value 覆盖。"
    elif method == "revised2009":
        if bun_mg_dl is None:
            raise ValueError("revised2009 需要提供 bun_mg_dl")
        denom = serum_creatinine_mgdl + 0.003 * bun_mg_dl - 0.024
        if denom <= 0:
            raise ValueError("revised2009 分母非正（Scr+0.003×BUN-0.024 必须 > 0）")
        egfr = _BEDSIDE_K * height_cm / denom
        formula = f"eGFR = 0.413×height/(Scr + 0.003×BUN − 0.024)"
        note = "含 BUN 修订 Schwartz 2009，对 eGFR<60 的儿童更准确。"
    else:  # bedside2009
        egfr = _BEDSIDE_K * height_cm / serum_creatinine_mgdl
        formula = "eGFR = 0.413×height/Scr"
        note = "床旁 Schwartz 2009（KDIGO 推荐默认式）。"

    egfr = round(egfr, 1)
    pediatric_caveat = (
        "注意：<2 岁婴儿 eGFR 参考范围低于年长儿（正常可低至 60–90），"
        "G1/G2 阈值在婴儿期需结合月龄与生长曲线解读。"
        if age_years < 2 else ""
    )
    return {
        "egfr": egfr,
        "unit": "ml/min/1.73m2",
        "method": method,
        "formula": formula,
        "pediatric_caveat": pediatric_caveat,
        "note": note,
    }


def classify_ckd(
    egfr: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
) -> dict:
    """KDIGO 2024 儿童 CKD 合并分期。

    身份来自部署注入的环境变量 A207_CALLER（P0-1：模型不可自证身份）。
    egfr 必需；白蛋白尿二选一（uacr 或 upcr）。返回 g、a、stage(GxAx)、description、risk_note。
    """
    enforce_read(MCP_NAME)
    if egfr < 0:
        raise ValueError("egfr 必须 >= 0")

    g = _egfr_to_g(egfr)
    a: Optional[str] = None
    albuminuria_source = ""
    if uacr_mg_g is not None:
        a = _acr_to_a(uacr_mg_g)
        albuminuria_source = "UACR"
    elif upcr_mg_g is not None:
        a = _pcr_to_a(upcr_mg_g)
        albuminuria_source = "UPCR"

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
        "stage": stage,
        "g": g,
        "a": a,
        "g_description": g_desc,
        "a_description": a_desc.get(a) if a else None,
        "albuminuria_source": albuminuria_source or None,
        "risk_note": risk_note,
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
        v = new_labs[metric]
        op = rule["operator"]
        if op == "gt":
            hit = v > rule["threshold"]
        elif op == "lt":
            hit = v < rule["threshold"]
        elif op == "between":
            low, high = rule["low"], rule["high"]
            hit = low <= v < high
        else:
            hit = False
        if not hit:
            return None
        return {
            "metric": metric,
            "observed": v,
            "threshold": (rule.get("low"), rule.get("high")) if op == "between"
                         else rule.get("threshold"),
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
                hit = -rule["high_pct"] < pct <= -rule["low_pct"]
            else:
                hit = pct <= -rule["threshold_pct"]
        if not hit:
            return None
        return {
            "metric": metric,
            "observed": round(pct, 1),
            "threshold": (f">= {rule['threshold_pct']}%" if direction == "up"
                          else f"<= -{rule['threshold_pct']}%"),
            "unit": "%",
        }
    return None


def evaluate_risk_rules(
    new_labs: Dict[str, float],
    prior_labs: Optional[Dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> Dict[str, Any]:
    """基于本轮新数据重评风险。

    入参：
      new_labs: 本轮 M2 拉回的新化验/测量（PCP 标准单位）。
                支持键：scr(mg/dL), k, p, hb(g/L), ca, ua(umol/L), egfr(ml/min/1.73m^2)
      prior_labs: 上轮同指标值（用于趋势类规则：R-01/R-08 肌酐、R-07 eGFR）。缺则趋势规则不触发。
      prior_level: 历史系统等级（L1/L2/L3）。**仅作对比输出，绝不兜底**。
      caller: 内部形参，缺省由部署注入的 A207_CALLER 解析（P0-1：模型不可自证身份）。

    返回：
      matched_rules: 命中规则明细列表
      overall_level: 最高等级（L1>L2>L3），无命中为 "none"
      prior_comparison: {prior_level, current_level, delta_note}
      level_correction_applied: True（声明范式已执行）
    """
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

    # 最高等级
    overall = "none"
    for m in matched:
        if _LEVEL_RANK[m["level"]] > _LEVEL_RANK[overall]:
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
        "matched_rules": matched,
        "overall_level": overall,
        "prior_comparison": {
            "prior_level": prior_level,
            "current_level": overall,
            "delta_note": delta_note,
        },
        "level_correction_applied": True,
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
                                  else f"up [{r['low_pct']}, {r['high_pct']})%")
        out.append(entry)
    return out


def explain_verdict(evaluation: Dict[str, Any]) -> List[Dict[str, Any]]:
    """把 evaluate_risk_rules 的结果翻成判定链路，供审计与医生复核。"""
    chain = []
    for m in evaluation["matched_rules"]:
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
    return chain

# ================================================================
# DAG: assess_clinical_status (v2.2/v2.3)
# ================================================================

def assess_clinical_status(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    new_labs: Optional[Dict[str, float]] = None,
    prior_labs: Optional[Dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> Dict[str, Any]:
    """一键评估 CKD 临床状态（eGFR + 分期 + 风险评分 DAG）。

    身份来自部署注入的环境变量 A207_CALLER（P0-1）。
    如果提供了 bun_mg_dl 则自动使用 revised2009 修订公式；若显式传入 k_value 则切换经典公式。
    """
    enforce_read(MCP_NAME)

    # 自动识别 Schwartz 方法
    method: EgfrMethod = "bedside2009"
    if bun_mg_dl is not None:
        method = "revised2009"
    elif k_value is not None:
        method = "classic"

    egfr_r = calc_egfr_schwartz(
        age_years=age_years,
        height_cm=height_cm,
        serum_creatinine_mgdl=serum_creatinine_mgdl,
        method=method,
        bun_mg_dl=bun_mg_dl,
        k_value=k_value,
    )
    ckd_r = classify_ckd(
        egfr=egfr_r["egfr"],
        uacr_mg_g=uacr_mg_g,
        upcr_mg_g=upcr_mg_g,
    )
    risk_r: Dict[str, Any] = {}
    if new_labs:
        labs = dict(new_labs)
        labs.setdefault("egfr", egfr_r["egfr"])
        risk_r = evaluate_risk_rules(
            new_labs=labs,
            prior_labs=prior_labs,
            prior_level=prior_level,
        )
    return {
        "ok": True,
        "egfr": egfr_r["egfr"],
        "egfr_unit": egfr_r["unit"],
        "egfr_method": egfr_r["method"],
        "ckd_stage": ckd_r["stage"],
        "ckd_g_stage": ckd_r["g"],
        "ckd_a_stage": ckd_r.get("a"),
        "risk_level": risk_r.get("overall_level", "none") if risk_r else "none",
        "risk_matched_rules": risk_r.get("matched_rules", []) if risk_r else [],
    }

# ---- M8: risk rules engine ----
