"""P4 评估计算域 MCP 服务：eGFR + CKD 分期 + 风险引擎 + DAG 一键评估。

合并自 M6 (a207-clinical-calc-mcp) + M8 (a207-risk-rules-mcp)。
v2.3 新增 DAG: assess_clinical_status（eGFR→分期→风险 一键判定）。
"""
from __future__ import annotations

import json
import logging

from typing import Any, Literal, Optional

from fastmcp import FastMCP

from a207_policy import CallerError

from .core import (
    assess_clinical_status,
    explain_verdict,
    list_rules,
)

mcp = FastMCP("CKDNutri-assessment-mcp")

# 2026-08-12（server 二次审查）：标准 logging 提升可观测性 + 异常分级归类——
# ① ValueError（core 层业务/参数校验，如 _require_finite/k_value/egfr 负数）归
# INVALID_INPUT 且 detail 保留（对编排层有明确语义）；② 未知系统异常（TypeError/
# KeyError/AttributeError/ZeroDivisionError 等内部 Code Bug）归 INTERNAL_ERROR 且
# detail **脱敏**（不裸暴露 str(exc) 内部实现，完整 StackTrace 仅留服务端日志）。
logger = logging.getLogger("CKDNutri-assessment-mcp")


def _invalid(exc: Exception) -> dict[str, Any]:
    if isinstance(exc, CallerError):
        # BUG-54（2026-08-12）：越权/身份未解析统一返回 FORBIDDEN 信封（与 care _guard /
        # clinical-data _guard_access 同格式），不再向上抛导致 500。此前本包 4 处裸调
        # enforce_read 的读工具越权即 500 崩溃。PermissionDenied 带 caller/action/reason，
        # CallerUnknown 缺字段时降级文案。
        # 2026-08-12（六审）：reason 三重保底（属性为空串时 getattr 默认值不生效）。
        # 2026-08-12（七审）：caller/action 亦做 or 保底（属性被显式置 None 时）。
        caller = getattr(exc, "caller", None) or "?"
        action = getattr(exc, "action", None) or "access"
        reason = getattr(exc, "reason", None) or str(exc) or "无明确原因"
        return {"ok": False, "error": "FORBIDDEN",
                "detail": f"caller={caller} 无权 {action}（{reason}）"}
    # BUG-52（2026-08-12）：内部数据错误归 INTERNAL_ERROR，避免误归 INVALID_INPUT
    # 2026-08-12（系统性审查，P1）：detail **脱敏**——FileNotFoundError/OSError 的 str
    # 含服务端绝对路径（如 /var/app/data/.../rules.json），原样返回泄露内部文件系统
    # 结构；完整异常（含路径）仅留服务端日志（与 care server 同口径）。
    if isinstance(exc, (FileNotFoundError, OSError, json.JSONDecodeError)):
        logger.warning("评估服务内部数据错误: %s", exc)
        return {"ok": False, "error": "INTERNAL_ERROR",
                "detail": "内部数据错误（error_code=ASSESS_DATA），详情见服务端日志"}
    if isinstance(exc, ValueError):
        # core 层业务/参数校验异常（有限性、范围、枚举等）——detail 对调用方有明确语义，保留
        return {"ok": False, "error": "INVALID_INPUT", "detail": str(exc)}
    # 未知系统异常 = 内部 Code Bug——归 INTERNAL_ERROR（编排层不应重试/误判入参问题），
    # detail 脱敏（不泄露内部实现），完整堆栈仅服务端日志。
    logger.error("评估服务未预期异常（内部 bug，error_code=ASSESS_UNKNOWN）", exc_info=exc)
    return {"ok": False, "error": "INTERNAL_ERROR",
            "detail": "评估服务内部错误（error_code=ASSESS_UNKNOWN），请查服务端日志"}


def main():
    mcp.run()


# ---- DAG: 一键评估（v2.2/v2.3 新增） ----

# 2026-08-12（server 二次审查）：显式 @mcp.tool(name=...) 锁定对外契约——FastMCP 默认用
# Python 函数名注册，函数重命名（如加/去 _tool 后缀）会静默改变 MCP Client 可见的工具名，
# 导致编排层 ToolNotFound。全项目（care/nutrition/content 等）工具名统一 _tool 后缀，
# 此处显式锁定现状契约名，防未来重命名破坏。

@mcp.tool(name="assess_clinical_status_tool")
def assess_clinical_status_tool(
    age_years: float,
    height_cm: float,
    serum_creatinine_mgdl: float,
    serum_creatinine_unit: Literal["mg_dL", "umol_L"] = "mg_dL",
    is_preterm: bool = False,
    sex: Optional[Literal["M", "F"]] = None,
    uacr_mg_g: Optional[float] = None,
    upcr_mg_g: Optional[float] = None,
    bun_mg_dl: Optional[float] = None,
    k_value: Optional[float] = None,
    method: Optional[Literal["bedside2009", "revised2009", "classic"]] = None,
    new_labs: Optional[dict[str, float]] = None,
    prior_labs: Optional[dict[str, float]] = None,
    prior_level: Optional[Literal["L1", "L2", "L3", "none"]] = None,
    dialysis_mode: Optional[Literal["none", "hemodialysis", "peritoneal"]] = None,
) -> dict[str, Any]:
    """一键评估：eGFR→CKD 分期→风险等级 完整链路。

    入参：患者年龄、身高、血清肌酐数值、可选蛋白尿/风险数据。
    serum_creatinine_mgdl（数值）的计量单位由 serum_creatinine_unit 指定：
    - mg_dL（默认）：数值即 mg/dL；
    - umol_L：数值为 µmol/L，服务端自动 ÷88.4 换算为 mg/dL 后计算（P1 get_labs
      返回 scr_umol_L，直接透传可避免 88 倍单位误差）。

    sex（M/F，可选）：经典 Schwartz ≥13y 青少年女性 k=0.55（男性 0.70），
    仅 method="classic" 时生效；不传或 method 非 classic 时不影响计算。
    is_preterm：<1y 早产儿在 method 未指定时自动采用经典式 k=0.33（防 eGFR 高估）。

    new_labs（本轮新化验，可选）支持两种键名（**务必使用以下标准 Key，勿自定义
    如 creatinine/potassium 等名称**）：
    - 规则短名：scr(mg/dL), k(mmol/L), p(mmol/L), hb(g/L), ca(mmol/L), na(mmol/L),
      ua(umol/L), egfr(ml/min/1.73m²), bun(mg/dL), uacr(mg/g), upcr(mg/g)；
    - P1 完整键名（自动换算）：scr_umol_L(÷88.4), k_mmol_L, p_mmol_L, hb_g_L,
      ca_mmol_L, na_mmol_L, ua_umol_L, egfr_ml_min, bun_mmol_L(×2.8),
      bun_mg_dl/bun_mg_dL, uacr_mg_g, upcr_mg_g（键名匹配不区分大小写）。
    规则引擎仅识别上述 Key；未识别的 Key 会被忽略导致对应规则静默跳过。
    prior_labs（上轮同指标历史值，可选）键名同 new_labs，用于动态趋势类规则
    （scr/egfr 环比），缺失则趋势规则不触发。

    出参：{egfr, egfr_unit, egfr_method, ckd_stage, risk_level, risk_matched_rules,
    risk_completeness, warnings} 全链路聚合。

    dialysis_mode（可选，A-S1 修复 2026-08-14）：none / hemodialysis / peritoneal——
    eGFR<15 且透析时分期标注 G5D（此前工具入口缺该参数，透析患儿经 MCP 评估
    永远返回 G5 而非 G5D，测试仅覆盖 core 层恰好掩盖）。
    """
    try:
        return assess_clinical_status(
            age_years, height_cm, serum_creatinine_mgdl,
            serum_creatinine_unit=serum_creatinine_unit,
            is_preterm=is_preterm,
            sex=sex,
            uacr_mg_g=uacr_mg_g, upcr_mg_g=upcr_mg_g,
            bun_mg_dl=bun_mg_dl, k_value=k_value,
            method=method,
            new_labs=new_labs, prior_labs=prior_labs, prior_level=prior_level,
            dialysis_mode=dialysis_mode,
        )
    except Exception as exc:
        return _invalid(exc)


# ---- M6/M8 底层拆分函数（BUG-27 修复，2026-08-12）----
# 需求 P4 只登记 3 个工具（assess_clinical_status_tool / list_rules_tool / explain_verdict_tool），
# calc_egfr_schwartz / classify_ckd / evaluate_risk_rules 是 DAG 的内部拆分步骤——
# 不再注册为 @mcp.tool（此前暴露给 LLM 增加 prompt token 消耗与误调用风险，LLM 可能绕开 DAG）。
# 编排层/测试直接 import CKDNutri_assessment_mcp.core 调用纯函数（2026-08-12：删除此前
# 带 _invalid 包装的 _tool 别名——无人引用，属死代码）。


@mcp.tool(name="list_rules_tool")
def list_rules_tool() -> dict[str, Any]:
    """列出所有活跃风险规则及其阈值（供临床复核规则引擎配置）。"""
    try:
        return list_rules()
    except Exception as exc:
        return _invalid(exc)


@mcp.tool(name="explain_verdict_tool")
def explain_verdict_tool(evaluation: dict[str, Any]) -> dict[str, Any]:
    """把评估结果翻成判定链路（哪条规则、什么数值、什么阈值），供审计与医生复核。

    evaluation 接受评估结果信封或裸 data 体：{ok: true, data: {matched_rules 或
    risk_matched_rules: [{id, name, level, observed, threshold, unit, description}],
    overall_level, prior_comparison}}。
    matched_rules 须为列表，每个元素须含 id/name/level/observed/threshold/unit/
    description 七个字段（缺字段或结构非法将报 INVALID_INPUT 而非静默降级）。
    evaluation 必须是 dict 对象（不接受 JSON 字符串等序列化形态）。
    """
    try:
        return explain_verdict(evaluation)
    except Exception as exc:
        return _invalid(exc)


if __name__ == "__main__":
    main()
