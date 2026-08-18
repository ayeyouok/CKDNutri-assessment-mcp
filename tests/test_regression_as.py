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
    # UPCR 路径（2026-08-18 修正，审查报告2）：UPCR 单独提供时按 KDIGO 蛋白尿分级
    # （儿科口径 A1<200/A2 200-500/A3>500 mg/g）映射 A 期，不再 a=None（此前 a=None
    # 迫使 LLM 套 UACR 30/300 表误判：UPCR 122 mg/g 被错判 A2，应为 A1）。
    r = core.classify_ckd(egfr=45.0, upcr_mg_g=500.0)
    assert r["data"]["a"] == "A2", r["data"]
    assert r["data"]["albuminuria_source"] == "UPCR", r["data"]
    # 用户场景复现：eGFR 35.1 + UPCR 122 mg/g → 应为 G3bA1（非 G3bA2）
    r2 = core.classify_ckd(egfr=35.1, upcr_mg_g=122.0)
    assert r2["data"]["stage"] == "G3bA1", r2["data"]
    assert r2["data"]["a"] == "A1", r2["data"]
    # mg/mmol 换算方向修正（P0-2）：442 mg/mmol × 8.84 ≈ 3907 mg/g → A3
    r3 = core.classify_ckd(egfr=35.1, upcr_mg_mmol=442.0)
    assert r3["data"]["a"] == "A3", r3["data"]
    assert r3["data"]["stage"] == "G3bA3", r3["data"]
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


def _tier_of(g, a):
    from CKDNutri_assessment_mcp import core
    note = core._risk_note(g, a)
    # 注意文案是"进展风险高/中等/低"（"风险高"非"高风险"），须按完整标记判断
    if "进展风险高" in note:
        return "高"
    if "进展风险中等" in note:
        return "中"
    return "低"


def test_risk_heatmap_kdigo_2024():
    """审查报告2 P0-1：_risk_note 须对齐 KDIGO 2024 进展风险热图（G×A）。

    原加性评分 g_rank+a_rank≥6 判"高"对非 G4+ 患者不可达，G1–G3b 高危格被错判。
    期望值（3 档：低/中/高）取自 canonical KDIGO 热图（绿/黄/橙-红）。
    """
    from CKDNutri_assessment_mcp import core
    expect = {
        ("G1", "A1"): "低", ("G1", "A2"): "中", ("G1", "A3"): "高",
        ("G2", "A1"): "低", ("G2", "A2"): "中", ("G2", "A3"): "高",
        ("G3a", "A1"): "低", ("G3a", "A2"): "中", ("G3a", "A3"): "高",
        ("G3b", "A1"): "中", ("G3b", "A2"): "高", ("G3b", "A3"): "高",
        ("G4", "A1"): "高", ("G4", "A2"): "高", ("G4", "A3"): "高",
        ("G5", "A1"): "高", ("G5", "A2"): "高", ("G5", "A3"): "高",
        ("G5D", "A1"): "高", ("G5D", "A2"): "高", ("G5D", "A3"): "高",
    }
    for (g, a), exp in expect.items():
        got = _tier_of(g, a)
        assert got == exp, f"{g}{a}: 期望 {exp}，实际 {got}"
    # a=None 不应捏造低风险，应提示补充 UACR
    assert "未计入" in core._risk_note("G3b", None), core._risk_note("G3b", None)


def test_normalize_labs_unit_conflict_rejected():
    """审查报告2 P1：UPCR/scr 双单位冲突须 fail-closed 拒绝（不再静默取先者）。"""
    from CKDNutri_assessment_mcp import core
    # scr 88.4× 冲突：scr_umol_L=176.8(→2.0) vs scr_mg_dl=0.5
    try:
        core._normalize_labs({"scr_umol_L": 176.8, "scr_mg_dl": 0.5})
    except ValueError:
        pass
    else:
        raise AssertionError("scr 双单位冲突应拒绝")
    # 一致值不误拒
    norm = core._normalize_labs({"scr_umol_L": 176.8, "scr_mg_dl": 2.0})
    assert abs(norm["scr"] - 2.0) < 1e-6, norm
    # P1 唯一 UPCR 契约键 upcr_mg_mmol 现在能被识别（不再静默丢失）
    norm2 = core._normalize_labs({"upcr_mg_mmol": 442.0})
    assert abs(norm2["upcr"] - 442.0 * 8.84) < 1e-3, norm2


def test_validate_rules_schema_nan_threshold():
    """审查报告2 P2：单边阈值 NaN/Inf 须被 schema 校验拦截（L1 规则不静默死亡）。"""
    from CKDNutri_assessment_mcp import core
    bad = {"rules": [{"id": "R-X", "name": "x", "level": "L1", "type": "absolute",
                      "metric": "k", "unit": "mmol/L",
                      "operator": "gt", "threshold": float("nan"),
                      "description": "d"}]}
    try:
        core._validate_rules_schema(bad)
    except RuntimeError:
        pass
    else:
        raise AssertionError("NaN 阈值应被 schema 校验拒绝")
    inf_bad = dict(bad)
    inf_bad["rules"] = [dict(bad["rules"][0], threshold=float("inf"))]
    try:
        core._validate_rules_schema(inf_bad)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Inf 阈值应被 schema 校验拒绝")


def test_require_finite_rejects_bool():
    """审查报告2 P2：_require_finite 须拒绝布尔（float(True)=1.0 会命中数值规则）。"""
    from CKDNutri_assessment_mcp import core
    for bad in (True, False):
        try:
            core._require_finite(bad, "k")
        except ValueError:
            pass
        else:
            raise AssertionError(f"布尔 {bad} 应被拒绝")


def test_prior_level_none_text():
    """审查报告2 P2：prior_level='none' 应标注'首次评估'，不读作'较历史等级 none 升高'。"""
    from CKDNutri_assessment_mcp import core
    r = core.evaluate_risk_rules(new_labs={"scr": 1.0}, prior_level="none")
    assert r["ok"] is True
    note = r["data"]["prior_comparison"]["delta_note"]
    assert "首次评估" in note, note
    assert "较历史等级 none" not in note, note


def test_explain_verdict_missing_data_envelope():
    """审查报告2 P2：信封校验对称——{ok:true} 缺 data 与 {ok:true,data:null} 同拒不静默。"""
    from CKDNutri_assessment_mcp import core
    for env in ({"ok": True}, {"ok": True, "data": None}):
        try:
            core.explain_verdict(env)
        except ValueError:
            pass
        else:
            raise AssertionError(f"异常信封 {env} 应被拒绝")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print(f"P4 A-S1/A-B1/A-B2/A-B5 REGRESSION OK（{len(fns)} 个用例）")




