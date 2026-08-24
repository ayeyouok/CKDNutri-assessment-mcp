"""十七审（2026-08-24）assessment-mcp 修复回归：A1 eGFR 上限 300 / A4 阈值文本统一。

覆盖：
- A1：calc_egfr_schwartz 上限 250→300——≥13y 男 classic k=0.70 高身材低肌酐
      （180cm/0.5）eGFR≈252 不再误拒；真实单位错配（>300）仍拒绝。
- A4：_eval_rule 阈值文本复用 _format_rule_threshold（单一事实源，含方向词），
      与 list_rules / explain_verdict 对齐。
"""
import os

os.environ.setdefault("A207_ENV", "test")
os.environ.setdefault("A207_ACCEPT_DEV_STORAGE", "1")
os.environ.setdefault("A207_CALLER", "doctor_assistant")
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
_POLICY = Path(__file__).resolve().parents[1].parents[1] / "a207-policy" / "src"
for p in (_SRC, _POLICY):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from CKDNutri_assessment_mcp import core  # noqa: E402


def test_egfr_high_normal_not_rejected():
    """A1：高身材低肌酐正常高滤过（≈252）不再误拒；>300 仍拒绝。"""
    res = core.calc_egfr_schwartz(age_years=15, height_cm=180,
                                  serum_creatinine_mgdl=0.5,
                                  method="classic", sex="M")
    egfr = res["data"]["egfr_raw"]
    assert 250 < egfr <= 300, egfr
    # 真实单位错配（scr 单位用错，eGFR 远超 300）仍拒绝
    try:
        core.calc_egfr_schwartz(age_years=15, height_cm=180,
                                serum_creatinine_mgdl=0.4,
                                method="classic", sex="M")
        raise AssertionError("eGFR>300 应拒绝")
    except ValueError:
        pass


def test_eval_rule_threshold_aligned_with_formatter():
    """A4：_eval_rule 阈值文本复用 _format_rule_threshold（含方向词）。"""
    down_rule = {"id": "R-x", "type": "trend_pct", "metric": "k",
                 "direction": "down", "threshold_pct": 30, "unit": "%"}
    r_down = core._eval_rule(down_rule, {"k": 1.0}, {"k": 2.0})  # -50% → 命中
    assert r_down is not None, "down 趋势应命中"
    assert r_down["threshold"] == core._format_rule_threshold(down_rule), r_down
    assert r_down["threshold"] == "down >= 30%", r_down

    up_rule = {"id": "R-y", "type": "trend_pct", "metric": "k",
               "direction": "up", "threshold_pct": 30, "unit": "%"}
    r_up = core._eval_rule(up_rule, {"k": 2.0}, {"k": 1.0})  # +100% → 命中
    assert r_up is not None, "up 趋势应命中"
    assert r_up["threshold"] == "up >= 30%", r_up
    assert r_up["threshold"] == core._format_rule_threshold(up_rule), r_up


if __name__ == "__main__":
    test_egfr_high_normal_not_rejected()
    test_eval_rule_threshold_aligned_with_formatter()
    print("ASSESSMENT 十七审 OK")
