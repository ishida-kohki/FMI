import sys
from importlib import import_module
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sample = import_module("sample")


def test_add() -> None:
    assert sample.add(2, 3) == 5


def test_sub() -> None:
    assert sample.sub(7, 4) == 3


def test_mul() -> None:
    assert sample.mul(6, 5) == 30


def test_div() -> None:
    assert sample.div(8, 2) == 4


def test_div_by_zero_raises() -> None:
    with pytest.raises(ValueError, match="division by zero"):
        sample.div(1, 0)
