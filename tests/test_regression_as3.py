"""十八审（2026-08-18）assessment-mcp 回归：审查报告 P1-1/P1-2/P1-3/P2-1/P2-2。

覆盖：
- P1-1：UPCR 来源 a_description 不再称"白蛋白尿"（蛋白尿交叉映射文案）
- P2-2：a_semantics 机器字段（albuminuria / proteinuria_crosswalk）
- P2-1：UPCR 双单位形参 canonicalize（同值合并放行 / 异值拒绝）
- P1-2：fully_evaluable 含 unknown_metrics 条件
- P1-3：assessment_status=FULL/PARTIAL 单字段
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

import pytest  # noqa: E402

from CKDNutri_assessment_mcp import core  # noqa: E402


def _classify(uacr=None, upcr=None, upcr_mmol=None, egfr=82.6):
    return core.classify_ckd(egfr=egfr, uacr_mg_g=uacr, upcr_mg_g=upcr,
                             upcr_mg_mmol=upcr_mmol)["data"]


def test_p11_upcr_a_description_proteinuria():
    """P1-1：仅 UPCR 来源时 a_description 为蛋白尿文案（非白蛋白尿）。"""
    d = _classify(upcr=300.0)  # A2（蛋白尿 200-500）
    assert d["a"] == "A2"
    assert d["a_semantics"] == "proteinuria_crosswalk"
    assert "蛋白尿" in d["a_description"] and "白蛋白尿" not in d["a_description"], d


def test_p11_uacr_a_description_albuminuria():
    """P1-1：UACR 来源时 a_description 保持白蛋白尿文案 + a_semantics=albuminuria。"""
    d = _classify(uacr=100.0)  # A2（白蛋白尿 30-300）
    assert d["a"] == "A2"
    assert d["a_semantics"] == "albuminuria"
    assert "白蛋白尿" in d["a_description"], d
    # 双指标时 UACR 优先（albuminuria_note 明示），语义仍为白蛋白尿
    d2 = _classify(uacr=100.0, upcr=800.0)
    assert d2["a_semantics"] == "albuminuria"
    assert d2["upcr_ignored"] is True


def test_p22_a_semantics_none_without_source():
    """P2-2：无白蛋白尿指标时 a=None、a_semantics=None（不捏造语义）。"""
    d = _classify()
    assert d["a"] is None
    assert d["a_semantics"] is None
    assert d["a_description"] is None


def test_p21_upcr_dual_unit_equivalent_merged():
    """P2-1：UPCR 双单位形参数值同值（100 mg/mmol == 884 mg/g）自动合并，不再误报。"""
    # DAG 入口（形参双单位同值 → canonicalize 合并，不再抛 M-5 双单位拒绝）
    r = core.assess_clinical_status(
        age_years=6, height_cm=120, serum_creatinine_mgdl=0.6,
        upcr_mg_g=884.0, upcr_mg_mmol=100.0)
    assert r["ok"] is True, r
    assert r["data"]["ckd_a_stage"] == "A3"  # 884 mg/g > 500 → A3
    # 直接 classify_ckd 双单位仍拒绝（直调防御不变）
    with pytest.raises(ValueError):
        core.classify_ckd(egfr=82.6, upcr_mg_g=884.0, upcr_mg_mmol=100.0)


def test_p21_upcr_dual_unit_conflict_rejected():
    """P2-1：UPCR 双单位形参数值异值显式拒绝（数据冲突，fail-closed）。"""
    with pytest.raises(ValueError):
        core.assess_clinical_status(
            age_years=6, height_cm=120, serum_creatinine_mgdl=0.6,
            upcr_mg_g=100.0, upcr_mg_mmol=100.0)  # 100 mg/g != 884 mg/g


def test_p12_fully_evaluable_includes_unknown_metrics():
    """P1-2：存在未识别指标（unknown_metrics 非空）时 fully_evaluable 必须为 False。"""
    r = core.assess_clinical_status(
        age_years=6, height_cm=120, serum_creatinine_mgdl=0.6,
        new_labs={"potassium": 5.5})  # "potassium" 不是规则短名 k
    d = r["data"]
    rc = d["risk_completeness"]
    assert rc["unknown_metrics"] == ["potassium"], rc
    assert rc["fully_evaluable"] is False, rc
    assert d["assessment_status"] == "PARTIAL", d


def test_p13_assessment_status_full():
    """P1-3：覆盖全部规则指标 + 趋势历史对照时 assessment_status=FULL。"""
    labs = {"ca": 2.4, "egfr": 82.6, "hb": 120, "k": 4.5, "na": 140,
            "p": 1.4, "scr": 0.6, "ua": 300}
    r = core.assess_clinical_status(
        age_years=6, height_cm=120, serum_creatinine_mgdl=0.6,
        new_labs=labs, prior_labs={"scr": 0.5, "egfr": 95.0})
    d = r["data"]
    assert d["risk_completeness"]["missing_metrics"] == [], d["risk_completeness"]
    assert d["risk_completeness"]["fully_evaluable"] is True
    assert d["assessment_status"] == "FULL", d
