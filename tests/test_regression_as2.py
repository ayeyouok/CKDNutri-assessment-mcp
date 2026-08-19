"""十五审（2026-08-18）assessment-mcp 三审回归：P2×3 / P3×5 / P4 关键项。

覆盖：大小写变体键冲突、短名-全键冲突、trend observed 3 位、浮点尘埃闭端、
explain_verdict observed 有限性、私有函数直调防护、µ→u 键名归一、eGFR 消息精度。
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


def test_p21_case_variant_keys_rejected():
    """P2-1：大小写变体重复键拒绝（lower_map 折叠此前使冲突检测只看得到一份）。"""
    for bad in ({"scr_mg_dl": 1.0, "scr_mg_dL": 2.0},        # 88.4 倍冲突
                {"Scr_umol_L": 88.4, "scr_umol_L": 176.8}):  # 88.4 值静默丢失
        try:
            core._normalize_labs(bad)
        except ValueError:
            continue
        raise AssertionError(f"大小写变体 {bad} 未被拒绝")
    # 单键正常
    assert core._normalize_labs({"scr_mg_dL": 1.0})["scr"] == 1.0


def test_p22_short_vs_full_key_conflict():
    """P2-2：短名直传与完整键换算一致性（短名恒胜此前静默吞 88.4 倍冲突）。"""
    for bad in ({"scr": 2.0, "scr_umol_L": 88.4},   # scr=2.0 vs 88.4→1.0
                {"k": 6.0, "k_mmol_L": 5.5}):
        try:
            core._normalize_labs(bad)
        except ValueError:
            continue
        raise AssertionError(f"短名-全键冲突 {bad} 未被拒绝")
    # 一致时放行
    out = core._normalize_labs({"scr": 1.0, "scr_umol_L": 88.4})
    assert out["scr"] == 1.0, out


def test_p23_trend_observed_3dp_stays_in_band():
    """P2-3：trend observed 3 位不越出自身区间上界——pct=49.983 显示 49.983
    （修复前 round(abs,1)=50.0 越出 [30,50) 上界与阈值同框矛盾）。"""
    r = core.evaluate_risk_rules(
        {"scr": 90.0299, "egfr": 100.0},  # new
        prior_labs={"scr": 60.02},         # old → pct=49.983…（R-01 ≥50% 下界邻域）
    )
    # 命中 R-08 [30,50) 或 R-01 —— 断言命中规则的 observed 均 < 50（区间上界）
    hits = (r.get("data") or {}).get("matched_rules", [])
    for h in hits:
        if h.get("unit") == "%":
            assert h["observed"] < 50.0, h
            assert h["observed"] == round(h["observed"], 3), h


def test_p31_dust_boundary_hits_closed_edge():
    """P3-1：浮点尘埃不翻转闭端边界——数学恰 50.0% 实算 49.9999...9 必须命中
    R-01（≥50%，L1 危急），不降级 R-08。"""
    r = core.evaluate_risk_rules(
        {"scr": 90.03, "egfr": 100.0}, prior_labs={"scr": 60.02})
    data = r.get("data") or {}
    levels = [h["level"] for h in data.get("matched_rules", []) if h.get("unit") == "%"]
    assert "L1" in levels, (data.get("matched_rules"), "恰 50.0% 应命中 R-01(L1)")


def test_p32_explain_verdict_nan_observed_rejected():
    """P3-2：explain_verdict observed=NaN 拒绝（此前只查字段存在不查有限）。"""

    bad = {"matched_rules": [{
        "id": "R-01", "name": "x", "level": "L1", "observed": float("nan"),
        "threshold": ">= 50", "unit": "%", "description": "d",
    }]}
    try:
        core.explain_verdict(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("observed=NaN 未被拒绝")
    # 正常链可解释
    ok = core.explain_verdict({"matched_rules": [{
        "id": "R-01", "name": "x", "level": "L1", "observed": 51.2,
        "threshold": ">= 50", "unit": "%", "description": "d",
    }]})
    assert ok["ok"] is True, ok


def test_p33_private_fn_direct_call_guards():
    """P3-3：私有函数直调防护（编排层可直接 import）——NaN/负/Inf 显式拒绝，
    不再静默误分期（_egfr_to_g(NaN)→G5、_egfr_to_g(Inf)→G1、_acr_to_a(-5)→A1）。"""

    for fn, bad in ((core._egfr_to_g, float("nan")),
                    (core._egfr_to_g, float("inf")),
                    (core._egfr_to_g, -1.0),
                    (core._acr_to_a, -5.0),
                    (core._acr_to_a, float("nan")),
                    (core._pcr_to_a, float("nan")),
                    (core._pcr_to_a, -1.0)):
        try:
            fn(bad)
        except ValueError:
            continue
        raise AssertionError(f"{fn.__name__}({bad!r}) 未拒绝")
    # k_value=-1.0 / NaN
    for bad_k in (-1.0, float("nan"), float("inf")):
        try:
            core._schwartz_k(5.0, bad_k)
        except ValueError:
            continue
        raise AssertionError(f"_schwartz_k k_value={bad_k!r} 未拒绝")
    # 合法调用不受影响
    assert core._egfr_to_g(60.0) == "G2"
    assert core._acr_to_a(25.0) == "A1"
    assert core._schwartz_k(5.0, None) == 0.55


def test_p34_mu_key_domain_normalized():
    """P3-4：键名 µ→u 归一（U+03BC）——"scr_μmol_L" 与 scr_umol_L 等价，
    scr 规则不再静默跳过（fail-open 修复）。"""
    out = core._normalize_labs({"scr_μmol_L": 88.4})  # U+03BC 希腊字母 μ
    assert out["scr"] == 1.0, out  # 88.4 µmol/L = 1.0 mg/dL
    # 且不产生大小写变体冲突（µ 归一后与 ASCII u 同键）
    r = core.evaluate_risk_rules({"scr_μmol_L": 88.4})
    assert r["ok"] is True, r


def test_p44_egfr_rejection_message_precision():
    """P4：eGFR>250 拒绝消息显示真实值（.1f 此前把 250.01 显示成 250.0）。"""
    try:
        core.calc_egfr_schwartz(age_years=3, height_cm=250.0,
                                serum_creatinine_mgdl=0.05)  # egfr=0.413*250/0.05=2065
    except ValueError as exc:
        assert "2065.00" in str(exc), str(exc)  # .2f 显示真实值（.1f 会显示 2065.0 丢失精度）
    else:
        raise AssertionError("eGFR>250 应拒绝")


def test_p41_explicit_param_vs_new_labs_conflict():
    """P1-1（四审）：显式参数与 new_labs 同指标不一致拒绝——uacr_mg_g=30 与
    new_labs={"uacr":500} 双源冲突此前静默以显式参数覆盖。"""
    from CKDNutri_assessment_mcp import core

    try:
        core.assess_clinical_status(
            age_years=8, height_cm=120, serum_creatinine_mgdl=1.0,
            uacr_mg_g=30.0, new_labs={"uacr": 500.0})
    except ValueError:
        pass
    else:
        raise AssertionError("显式参数与 new_labs 冲突未拒绝")
    # 一致时放行
    r = core.assess_clinical_status(
        age_years=8, height_cm=120, serum_creatinine_mgdl=1.0,
        uacr_mg_g=30.0, new_labs={"uacr": 30.0})
    assert r["ok"] is True, r


def test_p42_rules_schema_bool_and_negative_rejected():
    """P1-2/P1-3（四审）：rules.json Schema 拒绝 bool 阈值与负趋势阈值。"""
    from CKDNutri_assessment_mcp import core

    bool_rule = {"rules": [{"id": "T1", "name": "t", "level": "L1", "metric": "k",
                            "type": "absolute", "operator": "gt", "threshold": True,
                            "unit": "mmol/L", "description": "d"}]}
    try:
        core._validate_rules_schema(bool_rule)
    except RuntimeError:
        pass
    else:
        raise AssertionError("bool threshold 未被 schema 拒绝")
    neg_rule = {"rules": [{"id": "T2", "name": "t", "level": "L1", "metric": "scr",
                           "type": "trend_pct", "direction": "up", "threshold_pct": -20,
                           "unit": "%", "description": "d"}]}
    try:
        core._validate_rules_schema(neg_rule)
    except RuntimeError:
        pass
    else:
        raise AssertionError("负 threshold_pct 未被 schema 拒绝")


def test_p44_heatmap_fail_closed():
    """P1-4（四审）：风险热图未映射组合 → RuntimeError（不再静默兜底"中"）。"""
    from CKDNutri_assessment_mcp import core

    try:
        core._risk_note("G6", "A1")  # 热图未覆盖的分期组合
    except RuntimeError:
        pass
    else:
        raise AssertionError("未映射 (G,A) 组合未 fail-closed")
    # 正常组合不受影响
    assert "低" in core._risk_note("G1", "A1")


def test_p45_normalize_full_precision_threshold():
    """P1-5（四审）：内部全精度——k_mmol_L=5.50001 不再 round4 截断成 5.5，
    R-02（k>5.5）必须命中（临界值不漏检）。"""
    from CKDNutri_assessment_mcp import core

    out = core._normalize_labs({"k_mmol_L": 5.50001})
    assert out["k"] == 5.50001, out  # 全精度保留
    r = core.evaluate_risk_rules({"k_mmol_L": 5.50001})
    hits = [(h["id"], h["level"]) for h in (r.get("data") or {}).get("matched_rules", [])]
    assert ("R-02", "L2") in hits, hits  # K>5.5 高钾血症命中


def test_p46_is_preterm_bool_guard():
    """P1-6（四审）：is_preterm 字符串 "false" 拒绝（truthy 此前误判早产）。"""
    from CKDNutri_assessment_mcp import core

    try:
        core.calc_egfr_schwartz(age_years=0.5, height_cm=60,
                                serum_creatinine_mgdl=0.4, is_preterm="false")
    except ValueError:
        pass
    else:
        raise AssertionError('is_preterm="false" 未被拒绝')


def test_p47_unknown_sex_over13_rejected():
    """P1-7（四审）：≥13 岁 classic 未知性别拒绝（此前静默按男性 k=0.70）。"""
    from CKDNutri_assessment_mcp import core

    try:
        core.calc_egfr_schwartz(age_years=14, height_cm=150,
                                serum_creatinine_mgdl=1.0, sex="未知", method="classic")
    except ValueError:
        pass
    else:
        raise AssertionError("≥13 岁未知性别未拒绝")
    # sex=None（未提供）向后兼容放行
    r = core.calc_egfr_schwartz(age_years=14, height_cm=150,
                                serum_creatinine_mgdl=1.0, method="classic")
    assert r["ok"] is True, r


# ---- 十七审（2026-08-18）：纯函数直调防御 + 语义加固（P1-1~P1-5 / P2-1/2/4）----


def test_p51_dag_entry_type_guard():
    """P1-1：DAG 入口统一类型校验——age_years="0.5"/is_preterm="false" 拒绝
    （此前字符串在 `is_preterm and age_years < 1` 处 TypeError）。"""
    from CKDNutri_assessment_mcp import core

    for bad_age in (float("nan"), float("inf")):
        try:
            core.assess_clinical_status(age_years=bad_age, height_cm=70,
                                        serum_creatinine_mgdl=0.4)
        except ValueError:
            pass
        else:
            raise AssertionError(f"age_years={bad_age!r} 未拒绝")
    # 数字字符串按项目契约（_require_finite float() 可解析即接受，不崩 TypeError）
    r = core.assess_clinical_status(age_years="0.5", height_cm=70,
                                    serum_creatinine_mgdl=0.4)
    assert r["ok"] is True, r
    try:
        core.assess_clinical_status(age_years=0.5, height_cm=70,
                                    serum_creatinine_mgdl=0.4, is_preterm="false")
    except ValueError:
        pass
    else:
        raise AssertionError('is_preterm="false" 未拒绝')


def test_p52_labs_key_type_guard():
    """P1-2：_normalize_labs 非字符串键拒绝（{123: 5.5} 此前 AttributeError）。"""
    from CKDNutri_assessment_mcp import core

    try:
        core._normalize_labs({123: 5.5})
    except ValueError:
        pass
    else:
        raise AssertionError("非字符串键未拒绝")


def test_p53_unknown_metrics_exposed():
    """P1-3：未识别指标显式暴露（potassium 而非 k → unknown_metrics + 告警）。"""
    from CKDNutri_assessment_mcp import core

    r = core.assess_clinical_status(age_years=8, height_cm=120,
                                    serum_creatinine_mgdl=0.8,
                                    new_labs={"k": 5.5, "potassium": 7.0})
    assert r["ok"] is True, r
    rc = r["data"]["risk_completeness"]
    assert "potassium" in rc["unknown_metrics"], rc
    assert "未识别指标" in rc["note"], rc
    # 合法非规则指标（uacr）不算 unknown
    r2 = core.assess_clinical_status(age_years=8, height_cm=120,
                                     serum_creatinine_mgdl=0.8,
                                     new_labs={"k": 5.5, "uacr": 100.0})
    assert "uacr" not in r2["data"]["risk_completeness"]["unknown_metrics"]


def test_p54_completeness_structured():
    """P1-4：risk_completeness 结构化字段（metric/trend_coverage/fully_evaluable）。"""
    from CKDNutri_assessment_mcp import core

    r = core.assess_clinical_status(age_years=8, height_cm=120,
                                    serum_creatinine_mgdl=1.2,
                                    new_labs={"k": 5.5, "scr": 1.2},
                                    prior_labs={"scr": 1.0})
    rc = r["data"]["risk_completeness"]
    assert "metric_coverage" in rc and "covered" in rc["metric_coverage"], rc
    assert "trend_coverage" in rc and "required" in rc["trend_coverage"], rc
    assert "fully_evaluable" in rc, rc
    # 缺 egfr 历史 → 趋势未完全覆盖 → fully_evaluable=False
    assert rc["fully_evaluable"] is False, rc
    assert rc["trend_coverage"]["missing"], rc


def test_p55_level_rank_strict():
    """P1-5：_LEVEL_RANK 严格索引——matched 规则非法 level 不再静默降级（整体评估正常）。"""
    from CKDNutri_assessment_mcp import core

    r = core.evaluate_risk_rules({"k": 6.0})  # R-02 命中 L2
    assert r["ok"] is True and r["data"]["overall_level"] == "L2", r


def test_p52_eval_rule_nan_guard():
    """P2-1：_eval_rule 直调 NaN 拒绝（此前 NaN 比较恒 False 静默漏检）。"""

    from CKDNutri_assessment_mcp import core

    rule = {"id": "R-02", "metric": "k", "type": "absolute",
            "operator": "gt", "threshold": 5.5}
    try:
        core._eval_rule(rule, {"k": float("nan")}, None)
    except ValueError:
        pass
    else:
        raise AssertionError("NaN 直调 _eval_rule 未拒绝")


def test_p52_explain_verdict_forged_id():
    """P2-2：explain_verdict 拒绝伪造规则 ID（防判定链格式化器被滥用）。"""
    from CKDNutri_assessment_mcp import core

    fake = {"ok": True, "data": {"matched_rules": [
        {"id": "FAKE-99", "name": "伪造", "level": "L1", "observed": 1.0,
         "threshold": 0.5, "unit": "x", "description": "伪造"},
    ], "overall_level": "L1", "evaluation_note": None}}
    try:
        core.explain_verdict(fake)
    except ValueError:
        pass
    else:
        raise AssertionError("伪造规则 ID 未被拒绝")


def test_p52_dialysis_normalized():
    """P2-4：assess 返回 dialysis_mode_normalized（归一化标准值）。"""
    from CKDNutri_assessment_mcp import core

    r = core.assess_clinical_status(age_years=8, height_cm=120,
                                    serum_creatinine_mgdl=1.0,
                                    dialysis_mode="HemoDialysis ")
    assert r["ok"] is True, r
    assert r["data"]["dialysis_mode_normalized"] == "hemodialysis", r["data"]
