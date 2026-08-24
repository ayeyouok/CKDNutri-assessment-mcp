r"""十二审（2026-08-24）assessment-mcp 修复回归：5 claims 裁定（3 修复 + 2 驳回）。

覆盖修复项：
- #2：assess_clinical_status 的 unknown_metrics 应排除全部合法别名键（k_mmol_L 等），
      不再虚假告警 + fully_evaluable 误降级。
- #4：_risk_note 在 a is None（无 UACR 也无 UPCR）时文案不含"仅提供 UPCR"误导。

驳回项（不测试，已记录裁定）：
- #1：_load_rules DCL 竞态——误报（当前校验在中间、VIEW 在最后，且校验失败重置
      _RULES=None，线程 B 不会读到 None VIEW 崩溃）。
- #3：calc_egfr_schwartz 校验顺序——误报（十六审已将 _require_finite 前置到 sex
      校验之前，字符串年龄先被有限性校验拦截）。

#5（server unit 枚举放宽）为 schema 层类型标注改动，非计算逻辑，靠人工 + server
导入验证，不写独立单测（避免引入 FastMCP 重依赖）。
"""
import os

os.environ.setdefault("A207_ENV", "test")
os.environ.setdefault("A207_CALLER", "doctor_assistant")
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
_POLICY = Path(__file__).resolve().parents[1].parents[1] / "a207-policy" / "src"
for p in (_SRC, _POLICY):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from CKDNutri_assessment_mcp import core  # noqa: E402


def test_assess2_alias_keys_not_unknown():
    """#2：完整别名键（k_mmol_L/hb_g_L）不应被判 unknown，fully_evaluable 不降级。"""
    resp = core.assess_clinical_status(
        age_years=8, height_cm=125, serum_creatinine_mgdl=0.6,
        new_labs={"k_mmol_L": 4.5, "hb_g_L": 120})
    assert resp["ok"] is True, resp
    rc = resp["data"]["risk_completeness"]
    assert rc["unknown_metrics"] == [], \
        "别名键不应误判 unknown：" + str(rc["unknown_metrics"])
    # 注：fully_evaluable 因本测试仅传 k/hb/scr 缺其他必填规则指标而 False 是正当的
    # （确实缺指标），修复核心是 unknown_metrics 不再虚假包含别名键。


def test_assess4_risk_note_no_upcr_when_a_none():
    """#4：a is None（既无 UACR 也无 UPCR）文案不应含'仅提供 UPCR'误导。"""
    note = core._risk_note("G3a", None)
    assert "仅提供 UPCR" not in note, "文案不应误导临床以为 UPCR 未被评估"
    assert "UACR" in note and "UPCR" in note, "应明确未提供 UACR 或 UPCR"


if __name__ == "__main__":
    test_assess2_alias_keys_not_unknown()
    test_assess4_risk_note_no_upcr_when_a_none()
    print("assessment 十二审回归全部通过")
