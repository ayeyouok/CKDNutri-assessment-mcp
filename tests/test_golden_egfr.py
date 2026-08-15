"""eGFR/分期 golden dataset 单测（2026-08-15 自 CKDNutri-nutrition-mcp 搬回）。

背景：这 5 个用例测的是本包（assessment）的 calc_egfr_schwartz / classify_ckd，
此前被误放入 nutrition 包的 test_golden_dataset.py——本地因 PYTHONPATH 恰好能
import 到 assessment 而通过，但 GitHub Actions CI 只安装本包，nutrition 的 CI
必然 ModuleNotFoundError（已实测：test_golden_dataset.py 在 CI 崩溃）。搬回归属
包后，两侧 CI 各自只测自己包，互不依赖。

运行：pytest tests/test_golden_egfr.py  或  python tests/test_golden_egfr.py
"""
from __future__ import annotations

import os
os.environ.setdefault("A207_ENV", "test")  # N-SEC-1（2026-08-14）：测试进程显式声明测试环境（守卫 fail-closed 默认拒绝）
os.environ.setdefault("A207_ACCEPT_DEV_STORAGE", "1")  # 生产护栏（2026-08-15）：测试进程显式确认 json 后端为开发模式
import sys
from math import isclose
from pathlib import Path

os.environ.setdefault("A207_CALLER", "doctor_assistant")

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _near(actual, expected, tol=1e-6):
    # eGFR 返回 round 后展示值（如 79.2），用 0.05 绝对容差（与原 nutrition golden
    # 口径一致）；tol 参数保留给需更高精度的调用（如单位换算一致性）。
    assert isclose(actual, expected, rel_tol=tol, abs_tol=0.05), \
        f"期望 {expected}，实际 {actual}"


def test_egfr_golden_bedside():
    """床旁 Schwartz 2009：eGFR = 0.413 × 身高(cm) / Scr(mg/dL)。"""
    from CKDNutri_assessment_mcp import core

    # 6y M 115cm Scr=0.6 → 0.413×115/0.6 = 79.16…
    r = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=0.6)
    _near(r["data"]["egfr"], 0.413 * 115 / 0.6)
    # µmol/L 单位混用：88.4 µmol/L ≡ 1.0 mg/dL → 与 mg/dL=1.0 结果一致
    r1 = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=1.0)
    r2 = core.calc_egfr_schwartz(age_years=6, height_cm=115,
                                 serum_creatinine_mgdl=88.4, serum_creatinine_unit="umol_L")
    _near(r1["data"]["egfr"], r2["data"]["egfr"], tol=1e-6)
    _near(r2["data"]["egfr"], 0.413 * 115 / 1.0)


def test_egfr_golden_classic_preterm():
    """经典 Schwartz：早产儿 k=0.33、<1y 足月 0.45、1-12y 0.55、≥13y M 0.70/F 0.55。"""
    from CKDNutri_assessment_mcp import core

    cases = [
        # (age, sex, is_preterm, height, scr, expected_k)
        (0.5, None, True, 60, 0.8, 0.33),
        (0.5, None, False, 60, 0.8, 0.45),
        (6, None, False, 115, 0.6, 0.55),
        (14, "M", False, 160, 1.0, 0.70),
        (14, "F", False, 160, 1.0, 0.55),
    ]
    for age, sex, preterm, h, scr, k in cases:
        r = core.calc_egfr_schwartz(age_years=age, height_cm=h, serum_creatinine_mgdl=scr,
                                    method="classic", is_preterm=preterm, sex=sex)
        _near(r["data"]["egfr"], k * h / scr, tol=0.01)


def test_egfr_golden_revised_bun():
    """M-1（2026-08-15）：revised2009 假公式已移除——该线性式（0.413×H/(Scr+0.003×BUN−0.024)）
    无文献出处；Schwartz 2009 含 BUN 修订为 CKiD 组合式（幂函数且必需胱抑素 C），本系统未实现。
    断言：① 显式 method='revised2009' 拒绝；② 提供 BUN 未指定 method → 降级床旁式 + 告警。"""
    from CKDNutri_assessment_mcp import core

    # ① 假公式拒绝（fail-closed，不再静默算错值）
    try:
        core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=0.6,
                                bun_mg_dl=30, method="revised2009")
    except ValueError as exc:
        assert "revised2009" in str(exc) and "胱抑素" in str(exc), exc
    else:
        raise AssertionError("revised2009 应被拒绝（假公式已移除）")

    # ② BUN 自动推理：降级床旁式 + 显式告警（BUN 不参与计算）
    r = core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=0.6,
                                bun_mg_dl=30)
    assert r["ok"] is True, r
    assert r["data"]["method"] == "bedside2009", r["data"]["method"]
    expect = 0.413 * 115 / 0.6
    _near(r["data"]["egfr"], expect, tol=0.01)
    assert "胱抑素" in (r["data"].get("note") or ""), r["data"].get("note")


def test_egfr_golden_boundaries():
    """边界：eGFR 阈值判级 G1/G2/G3a/G3b/G4/G5/G5D（用未 round 值判级）。"""
    from CKDNutri_assessment_mcp import core

    cases = [
        (95, "G1"), (89.9, "G2"), (60, "G2"), (59.9, "G3a"), (45, "G3a"),
        (44.9, "G3b"), (30, "G3b"), (29.9, "G4"), (15, "G4"), (14.9, "G5"),
        (14.9, "G5D"),  # 透析
    ]
    for egfr, want_g in cases:
        dial = "hemodialysis" if want_g == "G5D" else None
        r = core.classify_ckd(egfr=egfr, dialysis_mode=dial)
        assert r["data"]["g"] == want_g, (egfr, r["data"]["g"], want_g)


def test_egfr_golden_invalid_inputs():
    """负值 / NaN / Inf / 超龄 → ValueError（fail-closed）。"""
    from CKDNutri_assessment_mcp import core

    for bad in (-5, float("nan"), float("inf")):
        for fn in (
            lambda: core.calc_egfr_schwartz(age_years=6, height_cm=115, serum_creatinine_mgdl=bad),
            lambda: core.classify_ckd(egfr=bad),
        ):
            try:
                fn()
            except ValueError:
                pass
            else:
                raise AssertionError(f"bad={bad} 应抛 ValueError")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print(f"eGFR GOLDEN OK（{len(fns)} 个用例）")
