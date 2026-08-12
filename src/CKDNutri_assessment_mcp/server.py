"""P4 评估计算域 MCP 服务：eGFR + CKD 分期 + 风险引擎 + DAG 一键评估。

合并自 M6 (a207-clinical-calc-mcp) + M8 (a207-risk-rules-mcp)。
v2.3 新增 DAG: assess_clinical_status（eGFR→分期→风险 一键判定）。
"""
from __future__ import annotations

import json

from typing import Any, Optional

from fastmcp import FastMCP

from a207_policy import CallerError

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
        # BUG-54（2026-08-12）：越权/身份未解析统一返回 FORBIDDEN 信封（与 care _guard /
        # clinical-data _guard_access 同格式），不再向上抛导致 500。此前本包 4 处裸调
        # enforce_read 的读工具越权即 500 崩溃。PermissionDenied 带 caller/action/reason，
        # CallerUnknown 缺字段时降级文案。
        return {"ok": False, "error": "FORBIDDEN",
                "detail": f"caller={getattr(exc, 'caller', '?')} 无权 {getattr(exc, 'action', 'access')}"
                          f"（{getattr(exc, 'reason', str(exc))}）"}
    # BUG-52（2026-08-12）：内部数据错误归 INTERNAL_ERROR，避免误归 INVALID_INPUT
    if isinstance(exc, (FileNotFoundError, OSError, json.JSONDecodeError)):
        return {"ok": False, "error": "INTERNAL_ERROR",
                "detail": f"内部数据错误：{exc}"}
    return {"ok": False, "error": "INVALID_INPUT", "detail": str(exc)}


def main():
    mcp.run()


# ---- DAG: 一键评估（v2.2/v2.3 新增） ----

@mcp.tool
def assess_clinical_status_tool(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    serum_creatinine_unit: str = "mg_dL",
    is_preterm: bool = False,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    method: Optional[str] = None,
    new_labs: Optional[dict[str, float]] = None,
    prior_labs: Optional[dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> dict[str, Any]:
    """一键评估：eGFR→CKD 分期→风险等级 完整链路。

    入参：患者年龄、身高、血清肌酐、可选蛋白尿/风险数据。
    serum_creatinine_unit: mg_dL（默认）或 umol_L——P1 get_labs 返回 scr_umol_L（µmol/L），
    传 umol_L 自动 ÷88.4 转 mg/dL（BUG-40 修复，防 88 倍误差）。
    出参：{egfr, ckd_stage, risk_level, risk_matched_rules} 全链路聚合。
    """
    try:
        return assess_clinical_status(
            age_years, height_cm, serum_creatinine_mgdl,
            serum_creatinine_unit=serum_creatinine_unit,
            is_preterm=is_preterm,
            uacr_mg_g=uacr_mg_g, upcr_mg_g=upcr_mg_g,
            bun_mg_dl=bun_mg_dl, k_value=k_value,
            method=method,
            new_labs=new_labs, prior_labs=prior_labs, prior_level=prior_level,
        )
    except Exception as exc:
        return _invalid(exc)


# ---- M6/M8 底层拆分函数（BUG-27 修复，2026-08-12）----
# 需求 P4 只登记 3 个工具（assess_clinical_status_tool / list_rules_tool / explain_verdict_tool），
# calc_egfr_schwartz / classify_ckd / evaluate_risk_rules 是 DAG 的内部拆分步骤——
# 不再注册为 @mcp.tool（此前暴露给 LLM 增加 prompt token 消耗与误调用风险，LLM 可能绕开 DAG）。
# 函数保留在模块内供编排层/测试直接 import 调用（经 CKDNutri_assessment_mcp.core）。

def calc_egfr_schwartz_tool(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    serum_creatinine_unit: str = "mg_dL",
    is_preterm: bool = False,
    method: str = "bedside2009",
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
) -> dict[str, Any]:
    """[内部] 估算肾小球滤过率（Schwartz 系列）。非 MCP 工具，供编排层 import 调用。
    serum_creatinine_unit: mg_dL / umol_L（BUG-40）；is_preterm: 早产儿 k=0.33（BUG-47）。"""
    try:
        return calc_egfr_schwartz(
            age_years, height_cm, serum_creatinine_mgdl,
            serum_creatinine_unit=serum_creatinine_unit,
            is_preterm=is_preterm,
            method=method, bun_mg_dl=bun_mg_dl, k_value=k_value,
        )
    except Exception as exc:
        return _invalid(exc)


def classify_ckd_tool(
    egfr: float,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
) -> dict[str, Any]:
    """[内部] KDIGO 2024 儿童 CKD 合并分期 GxA（含风险提示）。非 MCP 工具。"""
    try:
        return classify_ckd(egfr=egfr, uacr_mg_g=uacr_mg_g, upcr_mg_g=upcr_mg_g)
    except Exception as exc:
        return _invalid(exc)


def evaluate_risk_rules_tool(
    new_labs: dict[str, float],
    prior_labs: Optional[dict[str, float]] = None,
    prior_level: Optional[str] = None,
) -> dict[str, Any]:
    """[内部] 重评儿童 CKD 风险（L1/L2/L3 + 命中规则清单）。非 MCP 工具。"""
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
        return explain_verdict(evaluation)
    except Exception as exc:
        return _invalid(exc)


if __name__ == "__main__":
    main()
