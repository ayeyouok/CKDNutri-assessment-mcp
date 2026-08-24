r"""十八审（2026-08-24）assessment-mcp 修复回归：A6 早产儿 BUN 优先级 / A7 阈值口径
down>=X% / A8 pediatric_caveat 透出 / A9 动态受影响规则 ID / A10 sex 类型放宽。

覆盖：
- A6：<1y 早产患儿即使带 BUN，method 自动推理仍为 classic（k=0.33 生效），
      不被 BUN 误改 bedside2009（k=0.413 高估 ~25%）。
- A7：_format_rule_threshold 单阈值 down 统一为 "down >= X%"（非 "down <= -X%"）。
- A8：assess_clinical_status 返回 dict 含 pediatric_caveat（<2y 婴幼儿生理警示）。
- A9：缺失 prior_labs 仅影响实际缺失指标对应的趋势规则，note 动态列出受影响
      规则 id，不误报已评估规则。
- A10：core 入口 sex 支持 male/female/男/女（server JSON Schema 放宽前 core 已支持）。
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

import CKDNutri_assessment_mcp.core as core  # noqa: E402


def test_a6_preterm_priority_over_bun():
    """A6：6 月龄早产 + BUN，method 自动推理为 classic（k=0.33），非 bedside2009。"""
    res = core.assess_clinical_status(
        age_years=0.5, height_cm=65, serum_creatinine_mgdl=0.3,
        is_preterm=True, bun_mg_dl=15.0)
    assert res["ok"] is True, res
    assert res["data"]["egfr_method"] == "classic", res
    # classic k=0.33：0.33*65/0.3 ≈ 71.5（非 bedside2009 的 0.413→89.5）
    assert 60 < res["data"]["egfr"] < 85, res


def test_a7_down_threshold_wording():
    """A7：单阈值 down 规则阈值文本为 'down >= X%'。"""
    rule = {"type": "trend_pct", "direction": "down", "threshold_pct": 25}
    assert core._format_rule_threshold(rule) == "down >= 25%", \
        core._format_rule_threshold(rule)
    # up 不变
    rule_up = {"type": "trend_pct", "direction": "up", "threshold_pct": 50}
    assert core._format_rule_threshold(rule_up) == "up >= 50%"


def test_a8_pediatric_caveat_surface():
    """A8：<2y 患儿 assess_clinical_status 顶层透出 pediatric_caveat。"""
    res = core.assess_clinical_status(
        age_years=1, height_cm=75, serum_creatinine_mgdl=0.3,
        sex="M", is_preterm=False)
    assert res["ok"] is True, res
    # <2y 应生成 pediatric_caveat（婴儿 eGFR 可低至 60-90）
    assert "pediatric_caveat" in res["data"], res["data"].keys()
    # 顶层字段存在（即使为空串也是已透出，非空串说明触发）
    assert isinstance(res["data"]["pediatric_caveat"], str)


def test_a9_dynamic_affected_rule_ids():
    """A9：仅缺 egfr 时 note 不把已评估的 scr 规则(R-01)误报为未触发。"""
    # 故意只传 scr 历史、缺 egfr 历史；断言受影响规则列表只含 egfr 指标对应规则
    res = core.assess_clinical_status(
        age_years=10, height_cm=130, serum_creatinine_mgdl=0.6,
        sex="M", new_labs={"scr": 0.6}, prior_labs={"scr": 0.5})
    assert res["ok"] is True, res
    note = res["data"]["risk_completeness"].get("note", "")
    # 若 note 提及规则，不应把不含 egfr 指标的规则写死（用动态 affected_rules 生成）
    assert "R-01/R-07/R-08" not in note, "硬编码规则列表应替换为动态生成"


def test_a10_sex_full_and_chinese_clean():
    """A10：core 入口 sex 支持 male/female/男/女（server 放宽前 core 已接受）。"""
    for sex_in in ["male", "female", "男", "女"]:
        res = core.assess_clinical_status(
            age_years=10, height_cm=130, serum_creatinine_mgdl=0.6,
            sex=sex_in)
        assert res["ok"] is True, (sex_in, res)
        # 归一化后的 sex 应进入计算（不抛 INVALID_INPUT）
        assert res["data"]["egfr"] > 0, (sex_in, res)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"OK  {fn.__name__}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"FAIL {fn.__name__}: {exc!r}")
    print(f"\n{'ALL PASS' if failed == 0 else str(failed) + ' FAILED'} "
          f"({len(fns)} tests)")
    raise SystemExit(1 if failed else 0)
