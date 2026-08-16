# -*- coding: utf-8 -*-
"""A-S1/A-B1/A-B2/A-B5 回归测试（2026-08-14 修复后固化）。pytest + 直接运行双模式。"""
import os
os.environ.setdefault("A207_ENV", "test")  # N-SEC-1（2026-08-14）：测试进程显式声明测试环境（守卫 fail-closed 默认拒绝）
os.environ.setdefault("A207_ACCEPT_DEV_STORAGE", "1")  # 生产护栏（2026-08-15）：测试进程显式确认 json 后端为开发模式
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


def test_a_b6_classify_ckd_albuminuria():
    """F8（2026-08-16，第七轮审查）：classify_ckd 白蛋白尿 A 分期此前零测试
    （临床关键路径：A 分期决定随访频率）。KDIGO 2024：A1<30 / A2 30-300 / A3>300 mg/g。"""
    from CKDNutri_assessment_mcp import core

    def _a(uacr):
        r = core.classify_ckd(egfr=45.0, uacr_mg_g=uacr)
        assert r["ok"] is True, r
        return r["data"]["a"]

    assert _a(10) == "A1", _a(10)
    assert _a(30) == "A2", _a(30)    # KDIGO A2 = 30-300（含 300）
    assert _a(150) == "A2", _a(150)
    assert _a(300) == "A2", _a(300)
    assert _a(800) == "A3", _a(800)
    # UPCR 路径：KDIGO 2024 白蛋白尿 A 分期仅基于 UACR——UPCR 含球蛋白排泄与
    # 白蛋白不等价，设计上不映射 A 分期（a=None + 提示补充 UACR），断言该契约。
    r = core.classify_ckd(egfr=45.0, upcr_mg_g=500.0)
    assert r["data"]["a"] is None, r["data"]
    assert "UACR" in r["data"]["albuminuria_note"], r["data"]["albuminuria_note"]
    # 组合分期 GxAx
    r = core.classify_ckd(egfr=45.0, uacr_mg_g=150.0)
    assert r["data"]["stage"] == "G3aA2", r["data"]["stage"]

def test_a_b7_egfr_upper_bound_relaxed():
    """F1（2026-08-16，第七轮审查）：eGFR 上限 200→250——矮小+低肌酐正常变异
    （3 岁 100cm + scr 0.2 → eGFR≈206）不再被误拒；250 以上仍拒绝（录入错误）。"""
    from CKDNutri_assessment_mcp import core

    r = core.calc_egfr_schwartz(age_years=3, height_cm=100, serum_creatinine_mgdl=0.2)
    assert r["ok"] is True and r["data"]["egfr"] > 200, r
    try:
        core.calc_egfr_schwartz(age_years=10, height_cm=250, serum_creatinine_mgdl=0.05)
    except ValueError:
        pass
    else:
        raise AssertionError("eGFR>250 应拒绝（录入错误）")

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print(f"P4 A-S1/A-B1/A-B2/A-B5 REGRESSION OK（{len(fns)} 个用例）")




