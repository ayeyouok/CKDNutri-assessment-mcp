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


if __name__ == "__main__":
    test_server_imports()
    test_assess_and_explain()
    test_empty_labs_evaluation_note()
    test_rules_schema_validation()
    print("P4 SMOKE OK")
