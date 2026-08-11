"""P4 评估计算域 MCP 服务：eGFR + CKD 分期 + 风险引擎 + DAG 一键评估。

合并自 M6 (a207-clinical-calc-mcp) + M8 (a207-risk-rules-mcp)。
v2.3 新增 DAG: assess_clinical_status（eGFR→分期→风险 一键判定）。
"""
from __future__ import annotations

from typing import Any, Optional

from fastmcp import FastMCP

from ._policy import CallerError

from .core import (
    assess_clinical_status,
    calc_egfr_schwartz,
    classify_ckd,
    evaluate_risk_rules,
    explain_verdict,
    list_rules,
)

mcp = FastMCP("CKDNutri-assessment-mcp")


def _invalid(exc: Exception) -> dict[str, Any]:
    if isinstance(exc, CallerError):
        raise
    return {"ok": False, "error": "INVALID_INPUT", "detail": str(exc)}


def main():
    mcp.run()


# ---- DAG: 一键评估（v2.2/v2.3 新增） ----

@mcp.tool
def assess_clinical_status_tool(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    new_labs: Optional[dict] = None,
    prior_labs: Optional[dict] = None,
    prior_level: Optional[str] = None,
) -> dict[str, Any]:
    """一键评估：eGFR→CKD 分期→风险等级 完整链路。

    入参：患者年龄、身高、血清肌酐、可选蛋白尿/风险数据。
    出参：{egfr, ckd_stage, risk_level, risk_matched_rules} 全链路聚合。
    """
    try:
        return assess_clinical_status(
            age_years, height_cm, serum_creatinine_mgdl,
            uacr_mg_g=uacr_mg_g, upcr_mg_g=upcr_mg_g,
            bun_mg_dl=bun_mg_dl, k_value=k_value,
            new_labs=new_labs, prior_labs=prior_labs, prior_level=prior_level,
        )
    except Exception as exc:
        return _invalid(exc)


# ---- M6: eGFR & CKD 分期（保留独立工具供单点查询） ----

@mcp.tool
def calc_egfr_schwartz_tool(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    method: str = "bedside2009",
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
) -> dict[str, Any]:
    """估算肾小球滤过率（Schwartz 系列），CKD 临床助手独占。"""
    try:
        return calc_egfr_schwartz(
            age_years, height_cm, serum_creatinine_mgdl,
            method=method, bun_mg_dl=bun_mg_dl, k_value=k_value,
        )
    except Exception as exc:
        return _invalid(exc)


@mcp.tool
def classify_ckd_tool(
    egfr: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
) -> dict[str, Any]:
    """KDIGO 2024 儿童 CKD 合并分期 GxA（含风险提示）。"""
    try:
        return classify_ckd(egfr=egfr, uacr_mg_g=uacr_mg_g, upcr_mg_g=upcr_mg_g)
    except Exception as exc:
        return _invalid(exc)


# ---- M8: 风险规则引擎 ----

@mcp.tool
def evaluate_risk_rules_tool(
    new_labs: dict,
    prior_labs: Optional[dict] = None,
    prior_level: Optional[str] = None,
) -> dict[str, Any]:
    """重评儿童 CKD 风险（L1/L2/L3 + 命中规则清单）。"""
    try:
        return evaluate_risk_rules(
            new_labs=new_labs, prior_labs=prior_labs, prior_level=prior_level,
        )
    except Exception as exc:
        return _invalid(exc)


@mcp.tool
def list_rules_tool() -> dict[str, Any]:
    """列出所有活跃风险规则及其阈值。"""
    try:
        return list_rules()
    except Exception as exc:
        return _invalid(exc)


@mcp.tool
def explain_verdict_tool(evaluation: dict) -> dict[str, Any]:
    """把 evaluate_risk_rules 的返回翻成判定链路（哪条规则、哪个数值、什么阈值），供审计与医生复核。"""
    try:
        return {"chain": explain_verdict(evaluation)}
    except Exception as exc:
        return _invalid(exc)


if __name__ == "__main__":
    main()
