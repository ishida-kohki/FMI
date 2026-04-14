import sys
from importlib import import_module
from pathlib import Path


def test_sample_package_importable() -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    sample = import_module("sample")
    assert sample is not None
