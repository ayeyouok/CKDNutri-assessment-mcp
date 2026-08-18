# -*- coding: utf-8 -*-
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
    import math

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
    import math

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
