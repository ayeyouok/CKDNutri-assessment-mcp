# -*- coding: utf-8 -*-
"""A-S1/A-B1/A-B2/A-B5 回归测试（2026-08-14 修复后固化）。pytest + 直接运行双模式。"""
import os
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
_POLICY = Path(__file__).resolve().parents[1].parents[1] / "a207-policy" / "src"
for p in (_SRC, _POLICY):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
os.environ.setdefault("A207_CALLER", "doctor_assistant")


def test_a_s1_g5d_via_tool():
    """A-S1：assess_clinical_status_tool 透传 dialysis_mode → 透析患儿经 MCP 工具得 G5D。"""
    from CKDNutri_assessment_mcp import server

    r = server.assess_clinical_status_tool(
        age_years=6, height_cm=115, serum_creatinine_mgdl=8.0,  # eGFR≈5.9 <15
        dialysis_mode="peritoneal")
    assert r.get("ok") is True, r
    assert r["data"]["ckd_stage"] == "G5D", r["data"]["ckd_stage"]
    r2 = server.assess_clinical_status_tool(
        age_years=6, height_cm=115, serum_creatinine_mgdl=8.0)
    assert r2["data"]["ckd_stage"] == "G5", r2["data"]["ckd_stage"]


def test_a_b1_rules_cache_not_poisoned():
    """A-B1：schema 校验失败后 _RULES 重置（模块不永久损坏，修复配置后可恢复）。"""
    from CKDNutri_assessment_mcp import core

    core._reset_rules_cache()
    orig = core._validate_rules_schema

    def boom(doc):
        raise ValueError("模拟 schema 校验失败")

    core._validate_rules_schema = boom
    try:
        try:
            core._load_rules()
        except ValueError:
            pass
        assert core._RULES is None, "校验失败后缓存未重置（模块被投毒）"
    finally:
        core._validate_rules_schema = orig
    core._reset_rules_cache()
    assert core._load_rules() is not None, "重置后应能重新加载"


def test_a_b2_negative_labs_rejected():
    """A-B2：规则引擎物理区间校验（负值拒绝，与 eGFR 计算路径同深度）。"""
    from CKDNutri_assessment_mcp import core

    for bad in ({"hb": -10.0}, {"egfr": -10.0}, {"k": -1.0}):
        try:
            core.evaluate_risk_rules(new_labs=bad)
        except ValueError:
            continue
        raise AssertionError(f"负值 labs {bad} 应被拒绝")


def test_a_b5_short_key_case_insensitive():
    """A-B5：短名键大小写容错（HB=50 与 hb=50 等价，规则不静默跳过）。"""
    from CKDNutri_assessment_mcp import core

    r = core.evaluate_risk_rules(new_labs={"HB": 50.0})
    ids = {m["id"] for m in r["data"]["matched_rules"]}
    assert "R-04" in ids, ids
    assert r["data"]["overall_level"] == "L1", r["data"]["overall_level"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print(f"P4 A-S1/A-B1/A-B2/A-B5 REGRESSION OK（{len(fns)} 个用例）")
