"""十六审（2026-08-24）assessment-mcp 修复回归：#4 KDIGO 热图 / #2 阈值格式化 /
#5 classic+≥13y+sex=None 强制报错。

覆盖：
- #4：_risk_note G3aA1=中度、G3bA1=高风险（对齐 KDIGO 2012 官方热图）。
- #2：explain_verdict 阈值从规则库真实生成——between(R-10)/区间趋势(R-08)/
      单阈值趋势(up R-01 / down R-07) 均不再为 None 且 direction 正确。
- #5：calc_egfr_schwartz classic + age>=13 + sex=None + k_value=None 抛 ValueError
      （女性患儿 eGFR 静默高估 27% 的 fail-closed）。
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


def test_heatmap_g3a1_moderate_g3b1_high():
    """#4：KDIGO 2012 官方热图——G3aA1=中度(黄)、G3bA1=高风险(橙)。

    _RISK_HEATMAP 是 _risk_note 内局部变量，此处通过风险文案档位间接验证：
    G3aA1 → "进展风险中等"、G3bA1 → "进展风险高"。
    """
    # 绿区不受影响
    assert "进展风险低" in core._risk_note("G1", "A1")
    assert "进展风险低" in core._risk_note("G2", "A1")
    # #4 修复点：G3aA1 应为中度（此前误归"低"）
    assert "进展风险中等" in core._risk_note("G3a", "A1"), core._risk_note("G3a", "A1")
    # #4 修复点：G3bA1 应为高（此前误归"中"）
    assert "进展风险高" in core._risk_note("G3b", "A1"), core._risk_note("G3b", "A1")


def test_explain_verdict_threshold_from_rule_library():
    """#2：explain_verdict 阈值从规则库真实生成（between/区间/单阈值均不丢）。"""
    # between 规则 R-10
    r_between = core.explain_verdict({"matched_rules": [{
        "id": "R-10", "name": "x", "level": "L1", "observed": 4.0,
        "threshold": "?", "unit": "-", "description": "d"}]})
    assert r_between["ok"] is True
    assert "[" in r_between["data"]["chain"][0]["threshold"], r_between

    # 区间趋势 R-08（up [30,50)）
    r_range = core.explain_verdict({"matched_rules": [{
        "id": "R-08", "name": "x", "level": "L1", "observed": 40.0,
        "threshold": "?", "unit": "%", "description": "d"}]})
    assert "30%" in r_range["data"]["chain"][0]["threshold"], r_range

    # 单阈值 up R-01（up >=50%）
    r_up = core.explain_verdict({"matched_rules": [{
        "id": "R-01", "name": "x", "level": "L1", "observed": 51.2,
        "threshold": "?", "unit": "%", "description": "d"}]})
    assert r_up["data"]["chain"][0]["threshold"] == "up >= 50%", r_up

    # 单阈值 down R-07（down >= 25%）
    r_down = core.explain_verdict({"matched_rules": [{
        "id": "R-07", "name": "x", "level": "L2", "observed": 30.0,
        "threshold": "?", "unit": "%", "description": "d"}]})
    assert r_down["data"]["chain"][0]["threshold"] == "down >= 25%", r_down


def test_schwartz_classic_no_sex_female_overrate_rejected():
    """#5：classic + age>=13 + sex=None + k_value=None 抛 ValueError（防女性 eGFR 高估）。"""
    try:
        core.calc_egfr_schwartz(
            age_years=14, height_cm=160.0, serum_creatinine_mgdl=0.6,
            method="classic", sex=None, k_value=None)
    except ValueError as exc:
        assert "性别" in str(exc) or "sex" in str(exc).lower(), str(exc)
        return
    raise AssertionError("classic+≥13y+sex=None 未拒绝（女性 eGFR 将静默高估 27%）")


def test_schwartz_classic_explicit_k_value_ok_without_sex():
    """#5 对照：显式传 k_value 时即便 sex=None 也不应因缺性别拒绝（k 已由调用方决定）。"""
    res = core.calc_egfr_schwartz(
        age_years=14, height_cm=160.0, serum_creatinine_mgdl=0.6,
        method="classic", sex=None, k_value=0.55)
    assert res["ok"] is True, res
    # k=0.55 女性系数：eGFR = 0.55*160/0.6 ≈ 146.7
    assert 140 < res["data"]["egfr"] < 155, res
