"""二十审（2026-08-19）assessment-mcp 回归：P1-1/P1-2 + P2 建议项修复固化。

覆盖：
- P1-1：explain_verdict 拒绝伪造规则 metadata（真实 rule_id + 伪造 name/level/
  threshold/unit/description → 输出以规则库真实值覆盖）
- P1-2：_normalize_labs 单位换算后溢出（upcr_mg_mmol=1e308 → inf）拒绝
- P2-1：bedside2009 + k_value 显式拒绝（不静默忽略）
- P2-2：eGFR 分期边界固定回归（89.96 → 展示 90.0 但分期 G2；90.00 → G1 等）
- P2-3：G5D 透析优先（任意 eGFR + hemodialysis/peritoneal → G5D）
- P2-4：UPCR mg/mmol → mg/g（30 → 265.2 → A2）

pytest + 直接运行双模式（CI 逐文件 `python tests/test_*.py`，不依赖 pytest）。
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


def _expect_raises(exc_type, fn):
    try:
        fn()
    except exc_type:
        return
    raise AssertionError(f"expected {exc_type.__name__} to be raised")


# ---- P1-1：explain_verdict 防伪造 ----

def test_explain_verdict_rejects_forged_rule_metadata():
    """真实 rule_id + 伪造 metadata：输出必须用规则库真实值覆盖（防伪造审计）。"""
    # 取一个真实规则 id（规则库第一条第 R-01 类）
    rules = core._load_rules()["rules"]
    real_id = rules[0]["id"]
    real_rule = next(r for r in rules if r["id"] == real_id)
    fake_eval = {
        "ok": True,
        "data": {
            "matched_rules": [{
                "id": real_id,
                "name": "伪造规则名",
                "level": "L3",
                "observed": 1.0,
                "threshold": "伪造阈值",
                "unit": "伪造单位",
                "description": "伪造描述",
            }]
        }
    }
    r = core.explain_verdict(fake_eval)
    assert r["ok"] is True, r
    item = r["data"]["chain"][0]
    assert item["rule_name"] == real_rule["name"], (item, real_rule)  # 真实名
    assert item["level"] == real_rule["level"], (item, real_rule)     # 真实等级
    assert item["why"] == real_rule["description"], (item, real_rule)  # 真实描述
    assert item["unit"] == real_rule["unit"], (item, real_rule)       # 真实单位
    assert "伪造" not in str(item), item
    # observed（实测值）保留
    assert item["observed"] == 1.0, item


def test_explain_verdict_unknown_id_rejected():
    """未知规则 id 拒绝（既有行为回归）。"""
    bad = {"ok": True, "data": {"matched_rules": [
        {"id": "R-999", "name": "x", "level": "L1", "observed": 1.0,
         "threshold": 1, "unit": "u", "description": "d"}]}}
    _expect_raises(ValueError, lambda: core.explain_verdict(bad))


# ---- P1-2：单位换算后溢出 ----

def test_normalize_labs_rejects_conversion_overflow():
    """换算后 inf 必须拒绝（不进入规则计算）——放大因子（×8.84）溢出场景。

    scr_umol_L 是缩小因子（÷88.4，1e308 不溢出）故不在此断言；upcr_mg_mmol 的
    ×8.84 放大会把 1e308 变 inf（旧实现 inf 继续参与规则比较，isclose 语义异常）。
    """
    _expect_raises(ValueError, lambda: core._normalize_labs({"upcr_mg_mmol": 1e308}))
    # 正常值不受影响（upcr_mg_mmol → upcr ×8.84）
    out = core._normalize_labs({"upcr_mg_mmol": 30.0})
    assert abs(out["upcr"] - 265.2) < 1e-6, out


# ---- P2-1：bedside2009 + k_value ----

def test_bedside2009_k_value_rejected():
    """bedside2009 + k_value 显式拒绝（不静默忽略）。"""
    _expect_raises(ValueError, lambda: core.calc_egfr_schwartz(
        age_years=6, height_cm=120, serum_creatinine_mgdl=0.6,
        method="bedside2009", k_value=0.55))
    # classic + k_value 正常
    r = core.calc_egfr_schwartz(age_years=6, height_cm=120,
                                serum_creatinine_mgdl=0.6,
                                method="classic", k_value=0.55)
    assert r["ok"] is True, r


# ---- P2-2：eGFR 分期边界 ----

def test_egfr_stage_boundaries():
    """边界回归：分期用 egfr_raw（egfr 参数即原始值），四舍五入展示不得影响分期。"""
    cases = [(89.96, "G2"), (90.00, "G1"), (59.99, "G3a"), (60.00, "G2"),
             (44.99, "G3b"), (45.00, "G3a"), (29.99, "G4"), (30.00, "G3b"),
             (14.99, "G5"), (15.00, "G4")]
    for raw, want in cases:
        d = core.classify_ckd(egfr=raw)["data"]
        assert d["stage"] == want, (raw, d["stage"], want)
    # 89.96：raw 89.96 < 90 → 必须 G2（不得因展示四舍五入 90.0 误判 G1）；
    # classify_ckd 直接用 raw 分期，展示舍入在 eGFR 计算侧——此处验证分期口径。
    assert core.classify_ckd(egfr=89.96)["data"]["stage"] == "G2"


# ---- P2-3：G5D 透析优先 ----

def test_g5d_dialysis_priority():
    """透析状态必须始终优先 G5D（任意 eGFR 都不降级）。"""
    for egfr in (40.0, 100.0, 10.0):
        for dm in ("hemodialysis", "peritoneal"):
            d = core.classify_ckd(egfr=egfr, dialysis_mode=dm)["data"]
            assert d["stage"] == "G5D", (egfr, dm, d["stage"])


# ---- P2-4：UPCR 换算 ----

def test_upcr_unit_conversion():
    """UPCR 30 mg/mmol → 265.2 mg/g → A2（交叉映射）。"""
    # 换算本身在 _normalize_labs（×8.84 到 upcr_mg_g）
    out = core._normalize_labs({"upcr_mg_mmol": 30.0})
    assert abs(out["upcr"] - 265.2) < 1e-6, out
    # classify 分期（A2：200-500 mg/g）+ 语义标注
    d = core.classify_ckd(egfr=82.6, upcr_mg_g=265.2)["data"]
    assert d["a"] == "A2", d
    assert d["a_semantics"] == "proteinuria_crosswalk", d
    # 双单位同时提供 → 拒绝（DAG 入口同值合并放行，classify 直调拒绝）
    _expect_raises(ValueError, lambda: core.classify_ckd(
        egfr=82.6, upcr_mg_g=265.2, upcr_mg_mmol=30.0))


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print(f"AS4 REGRESSION OK（{len(fns)} 个用例）")
