"""Reusable scientific utilities for this repository."""


def add(a: float, b: float) -> float:
    return a + b


def sub(a: float, b: float) -> float:
    return a - b


def mul(a: float, b: float) -> float:
    return a * b


def div(a: float, b: float) -> float:
    if b == 0:
        msg = "division by zero"
        raise ValueError(msg)
    return a / b


__all__ = ["add", "div", "mul", "sub"]
