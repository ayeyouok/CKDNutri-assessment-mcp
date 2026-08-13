"""P4 冒烟自测：导入 server 不报错 + DAG 评估与判定解释可调用。

运行：pytest tests/test_import_smoke.py  (或 python tests/test_import_smoke.py)
依赖：a207-policy 已随 pip install -e . 安装。
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("A207_CALLER", "doctor_assistant")

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def test_server_imports():
    """导入 server 不可抛错（回归：explain_verdict_tool 曾错传关键字参数致 TypeError）。"""
    mod = importlib.import_module("CKDNutri_assessment_mcp.server")
    assert mod.mcp is not None


def test_assess_and_explain():
    from CKDNutri_assessment_mcp import core

    dag = core.assess_clinical_status(
        age_years=6, height_cm=115, serum_creatinine_mgdl=0.6,
        uacr_mg_g=40, new_labs={"scr": 0.6, "k": 4.5, "hb": 105}, prior_level="L2",
    )
    # BUG-15：DAG 返回统一 {ok, data} 信封；BUG-21：透出 prior_comparison
    assert dag.get("ok") is True
    assert "ckd_stage" in dag["data"] and "risk_level" in dag["data"]
    assert "prior_comparison" in dag["data"], "DAG 应透出 prior_comparison（BUG-21）"
    # 四审：DAG 透出计算层警示（egfr_note / ckd_risk_note）
    assert "egfr_note" in dag["data"] and dag["data"]["egfr_note"], "DAG 应透出 egfr_note"
    assert "ckd_risk_note" in dag["data"], "DAG 应透出 ckd_risk_note"

    ev = core.evaluate_risk_rules(new_labs={"scr": 0.6, "hb": 90},
                                  prior_labs={"scr": 0.55, "hb": 105}, prior_level="L2")
    chain = core.explain_verdict(ev)
    assert chain.get("ok") is True and len(chain["data"]["chain"]) >= 1
    # 四审：正常输入不触发空输入提示
    assert ev["data"].get("evaluation_note") is None, ev["data"]


def test_empty_labs_evaluation_note():
    """四审（2026-08-12）回归：evaluate_risk_rules 空输入显式提示（防"无风险"假象）。"""
    from CKDNutri_assessment_mcp import core

    ev = core.evaluate_risk_rules(new_labs={})
    assert ev.get("ok") is True and ev["data"]["overall_level"] == "none"
    assert "未提供任何化验指标" in (ev["data"].get("evaluation_note") or ""), ev["data"]


def test_rules_schema_validation():
    """四审（2026-08-12）回归：规则库结构校验 fail-closed——非法配置加载期拒绝。"""
    from CKDNutri_assessment_mcp import core

    def _rejects(rules_doc, label):
        try:
            core._validate_rules_schema(rules_doc)
        except ValueError:
            return
        raise AssertionError(f"期望 {label} 抛 ValueError")

    _rejects({"rules": []}, "空 rules")
    _rejects({"rules": [{"id": "X"}]}, "缺必填键")
    _rejects({"rules": [{"id": "X", "name": "n", "level": "L9", "type": "absolute",
                         "metric": "k", "unit": "mmol/L", "description": "d",
                         "operator": "gt", "threshold": 5.5}]}, "level 非法")
    _rejects({"rules": [{"id": "X", "name": "n", "level": "L1", "type": "absolute",
                         "metric": "k", "unit": "mmol/L", "description": "d",
                         "operator": "sideways", "threshold": 5.5}]}, "operator 非法")
    _rejects({"rules": [{"id": "X", "name": "n", "level": "L1", "type": "absolute",
                         "metric": "k", "unit": "mmol/L", "description": "d",
                         "operator": "between", "low": 5.0, "high": 5.0}]}, "between 边界")
    _rejects({"rules": [{"id": "X", "name": "n", "level": "L1", "type": "trend_pct",
                         "metric": "scr", "unit": "%", "description": "d",
                         "direction": "sideways", "threshold_pct": 50}]}, "direction 非法")
    _rejects({"rules": [{"id": "X", "name": "n", "level": "L1", "type": "unknown_type",
                         "metric": "scr", "unit": "%", "description": "d"}]}, "type 非法")
    # 合法规则不误拦
    core._validate_rules_schema({"rules": [
        {"id": "R-01", "name": "n", "level": "L1", "type": "absolute", "metric": "k",
         "unit": "mmol/L", "description": "d", "operator": "gt", "threshold": 5.5},
        {"id": "R-02", "name": "n", "level": "L2", "type": "trend_pct", "metric": "scr",
         "unit": "%", "description": "d", "direction": "up", "low_pct": 30, "high_pct": 50},
    ]})
    # 现有内置 rules.json 必须通过自身校验（防配置漂移回归）
    rules_doc = core._load_rules()
    core._validate_rules_schema(dict(rules_doc))


def test_s4_unauthorized_nan_unit():
    """S4（2026-08-13）补全：越权 / NaN / 单位换算自动化用例。"""
    from math import isclose

    from a207_policy import PermissionDenied

    from CKDNutri_assessment_mcp import core

    # ① 越权：parent_assistant 对 P4 矩阵 = ACCESS_NONE，直接调计算工具必须被拒
    os.environ["A207_CALLER"] = "parent_assistant"
    try:
        try:
            core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=0.6)
        except PermissionDenied:
            pass
        else:
            raise AssertionError("家长调用 calc_egfr_schwartz 应抛 PermissionDenied")
    finally:
        os.environ["A207_CALLER"] = "doctor_assistant"

    # ② NaN/Inf：scr 非有限值拒绝（六审已修 _normalize_scr 先 _require_finite）
    for bad in (float("nan"), float("inf")):
        try:
            core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"scr={bad} 应抛 ValueError")
    try:
        core.scr_umol_to_mgdl(float("nan"))
    except ValueError:
        pass
    else:
        raise AssertionError("scr_umol_to_mgdl(NaN) 应抛 ValueError")

    # ③ 单位换算：88.4 µmol/L ≡ 1.0 mg/dL（eGFR 结果一致）；显式换算函数正确
    a = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=1.0)
    b = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=88.4,
                                serum_creatinine_unit="umol_L")
    assert a["ok"] is True and b["ok"] is True
    assert isclose(a["data"]["egfr"], b["data"]["egfr"], rel_tol=1e-9), (a, b)
    assert isclose(core.scr_umol_to_mgdl(88.4), 1.0, rel_tol=1e-9)


def test_plausibility_range_and_explain_null():
    """2026-08-13（assessment 专项审查）回归：
    ① 临床合理性范围——age/height/scr 超物理区间拒绝（此前仅有限性/非负校验）；
    ② explain_verdict 对 {ok:true, data:null} 异常信封显式报错（此前静默当"无规则命中"）。"""
    from CKDNutri_assessment_mcp import core

    def _raises(fn, label):
        try:
            fn()
        except ValueError:
            return
        raise AssertionError(f"期望 {label} 抛 ValueError")

    # ① 合理性范围
    _raises(lambda: core.calc_egfr_schwartz(age_years=200, height_cm=115,
                                            serum_creatinine_mgdl=0.6), "age=200 超上限")
    _raises(lambda: core.calc_egfr_schwartz(age_years=6, height_cm=9999,
                                            serum_creatinine_mgdl=0.6), "height=9999 超上限")
    _raises(lambda: core.calc_egfr_schwartz(age_years=6, height_cm=115,
                                            serum_creatinine_mgdl=1e-6), "scr=1e-6 超下限")
    # 合法边界不误拦
    r = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=0.6)
    assert r["ok"] is True

    # ② explain_verdict 异常信封 fail-closed
    _raises(lambda: core.explain_verdict({"ok": True, "data": None}),
            "{ok:true, data:null} 信封")
    # 正常信封与裸 data 仍兼容
    ev = core.evaluate_risk_rules(new_labs={"scr": 0.6, "hb": 90})
    assert core.explain_verdict(ev)["ok"] is True
    assert core.explain_verdict(ev["data"])["ok"] is True


def test_s1_dialysis_g5d():
    """S1 修复（2026-08-13）回归：透析患儿 eGFR<15 → G5D（care 层据此给每月随访）。"""
    from CKDNutri_assessment_mcp import core

    # 透析 + eGFR<15 → G5D（hemodialysis / peritoneal / 大小写容错）
    for mode in ("hemodialysis", "peritoneal", "Hemodialysis"):
        r = core.classify_ckd(egfr=12.0, dialysis_mode=mode)
        assert r["data"]["g"] == "G5D", (mode, r)
    # 非透析 eGFR<15 → G5（G5D 不误报）
    assert core.classify_ckd(egfr=12.0)["data"]["g"] == "G5"
    assert core.classify_ckd(egfr=12.0, dialysis_mode="none")["data"]["g"] == "G5"
    # 透析但 eGFR≥15 → 仍按数值分期（不透支 G5D）
    assert core.classify_ckd(egfr=18.0, dialysis_mode="hemodialysis")["data"]["g"] == "G4"
    # dialysis_mode 非字符串 → ValueError（fail-closed）
    try:
        core.classify_ckd(egfr=12.0, dialysis_mode=123)
    except ValueError:
        pass
    else:
        raise AssertionError("dialysis_mode=123 应抛 ValueError")
    # DAG 端到端透传：透析 G5D / 非透析 G5
    dag = core.assess_clinical_status(age_years=6, height_cm=115,
                                      serum_creatinine_mgdl=4.0,
                                      dialysis_mode="hemodialysis")
    assert dag["data"]["ckd_g_stage"] == "G5D", dag["data"]
    dag2 = core.assess_clinical_status(age_years=6, height_cm=115,
                                       serum_creatinine_mgdl=4.0)
    assert dag2["data"]["ckd_g_stage"] == "G5", dag2["data"]


def test_contract_boundaries():
    """2026-08-13（assessment 二审 #7）契约测试：BUG-45 键名别名 / BUG-40 µmol 换算 /
    BUG-58 NaN / sex k / preterm / R-01/R-07/R-08 边界 / _invalid 分级 / explain 结构。

    此前这些修复全靠注释声明、无代码锁定——本轮实测全部行为正确后固化为回归。
    """
    from math import isclose

    from CKDNutri_assessment_mcp import core
    from CKDNutri_assessment_mcp.server import _invalid

    # ① BUG-45 键名别名：scr mg/dL 双变体 → R-01 命中（此前缺失，编排层直传静默漏报 AKI）
    # 数值用精确 +50%（1.5/1.0）——1.2/0.8 浮点误差 49.999…% 恰落阈值下，非别名问题。
    for key in ("scr_mg_dL", "scr_mg_dl"):
        r = core.evaluate_risk_rules(new_labs={key: 1.5}, prior_labs={key: 1.0})
        assert r["ok"] is True, r
        ids = {m["id"] for m in r["data"]["matched_rules"]}
        assert "R-01" in ids, f"{key} 别名未触发 R-01（BUG-45 回归）"

    # ② BUG-40 µmol/L 完整键名 → 等价命中（88.4 换算）
    r = core.evaluate_risk_rules(new_labs={"scr_umol_L": 88.4 * 1.5},
                                 prior_labs={"scr_umol_L": 88.4 * 1.0})
    assert "R-01" in {m["id"] for m in r["data"]["matched_rules"]}, "umol_L 未触发 R-01"

    # ③ BUG-58 NaN 全路径拒绝（evaluate 入口）
    try:
        core.evaluate_risk_rules(new_labs={"k": float("nan")})
    except ValueError:
        pass
    else:
        raise AssertionError("evaluate NaN 应抛 ValueError")

    # ④ sex k：classic ≥13y 女性 k=0.55 vs 男性 0.70（eGFR 比值锁定）
    f = core.calc_egfr_schwartz(14, 160, 1.0, method="classic", sex="F")["data"]["egfr"]
    m = core.calc_egfr_schwartz(14, 160, 1.0, method="classic", sex="M")["data"]["egfr"]
    assert isclose(f / m, 0.55 / 0.70, rel_tol=1e-3), (f, m)

    # ⑤ preterm：<1y is_preterm k=0.33 vs 足月 0.45
    p = core.calc_egfr_schwartz(0.5, 60, 0.8, method="classic", is_preterm=True)["data"]["egfr"]
    t = core.calc_egfr_schwartz(0.5, 60, 0.8, method="classic", is_preterm=False)["data"]["egfr"]
    assert isclose(p, 0.33 * 60 / 0.8, rel_tol=1e-2), p
    assert isclose(t, 0.45 * 60 / 0.8, rel_tol=1e-2), t

    # ⑥ R-01/R-08 边界：+50% 恰好命中 R-01 且不落 R-08；+30% 命中 R-08 且不落 R-01
    ids50 = {m["id"] for m in core.evaluate_risk_rules(
        new_labs={"scr": 1.5}, prior_labs={"scr": 1.0})["data"]["matched_rules"]}
    assert "R-01" in ids50 and "R-08" not in ids50, ids50
    ids30 = {m["id"] for m in core.evaluate_risk_rules(
        new_labs={"scr": 1.3}, prior_labs={"scr": 1.0})["data"]["matched_rules"]}
    assert "R-08" in ids30 and "R-01" not in ids30, ids30

    # ⑦ R-07 down 单阈值 -25%：恰好命中；-24.8% 不命中
    assert "R-07" in {m["id"] for m in core.evaluate_risk_rules(
        new_labs={"egfr": 45.0}, prior_labs={"egfr": 60.0})["data"]["matched_rules"]}
    assert "R-07" not in {m["id"] for m in core.evaluate_risk_rules(
        new_labs={"egfr": 45.1}, prior_labs={"egfr": 60.0})["data"]["matched_rules"]}

    # ⑧ _invalid 错误分级（BUG-52/54）：ValueError→INVALID_INPUT、数据错误→INTERNAL_ERROR+脱敏
    assert _invalid(ValueError("入参错"))["error"] == "INVALID_INPUT"
    assert _invalid(FileNotFoundError("/var/app/data/rules.json"))["error"] == "INTERNAL_ERROR"
    assert "/var/app" not in str(_invalid(FileNotFoundError("/var/app/data/rules.json")))

    # ⑨ explain_verdict 结构校验：残缺 matched_rules 显式报错（不静默当"无命中"）
    try:
        core.explain_verdict({"ok": True, "data": {"matched_rules": [{"id": "R-01"}]}})
    except ValueError:
        pass
    else:
        raise AssertionError("残缺 matched_rules 应抛 ValueError")


if __name__ == "__main__":
    test_server_imports()
    test_assess_and_explain()
    test_empty_labs_evaluation_note()
    test_rules_schema_validation()
    test_s4_unauthorized_nan_unit()
    test_plausibility_range_and_explain_null()
    test_s1_dialysis_g5d()
    test_contract_boundaries()
    print("P4 SMOKE OK")
