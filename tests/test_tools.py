"""CKDNutri-assessment-mcp（P4）自测入口（对齐 P1 test_tools.py 模式）。

N3 修复（2026-08-13）：assessment 此前只有 test_import_smoke.py 单一文件，
无统一测试运行器——与其他包（P1 test_tools.py / P2-P5 smoke）入口风格不一致，
CI 也难以统一调用。本文件把 smoke 全部用例聚合为一个可 import 的入口：

    python tests/test_tools.py            # 跑全部用例
    pytest tests/test_tools.py            # 或 pytest 模式

用例全部来自 test_import_smoke.py（含 S1 G5D / S4 越权 / 契约边界等），
保持单一事实源：本文件仅 re-export 并触发执行，不重复断言逻辑。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("A207_CALLER", "doctor_assistant")

import test_import_smoke as _smoke  # noqa: E402,F401  —— 执行全部 test_* 注册（pytest 模式）


def main() -> int:
    """python 模式：显式逐个执行 smoke 全部 test_* 函数。"""
    fns = sorted(
        (name, fn) for name, fn in vars(_smoke).items()
        if name.startswith("test_") and callable(fn)
    )
    for name, fn in fns:
        fn()
    print(f"P4 TOOLS OK（{len(fns)} 个用例）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
