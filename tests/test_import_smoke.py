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
    assert "ckd_stage" in dag and "risk_level" in dag

    ev = core.evaluate_risk_rules(new_labs={"scr": 0.6, "hb": 90},
                                  prior_labs={"scr": 0.55, "hb": 105}, prior_level="L2")
    chain = core.explain_verdict(ev)
    assert isinstance(chain, list) and len(chain) >= 1


if __name__ == "__main__":
    test_server_imports()
    test_assess_and_explain()
    print("P4 SMOKE OK")
