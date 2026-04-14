# AGENTS Guide

This repository is a Python-first scientific research skeleton for physics simulation and data analysis.

Use this file as guidance when contributing code, analysis workflows, and automation.

## Project Direction

- Prioritize numerical correctness, reproducibility, and clear experiment tracking.
- Keep reusable scientific logic in importable Python modules.
- Use interactive notebooks for exploration and presentation, not as the source of core algorithms.

## Recommended Repository Layout

- `src/sample/` for reusable package code.
- `notebooks/` for marimo apps (`*.py`).
- `tests/` for automated tests.
- `docs/` for project documentation.
- `work/` for large binary outputs, temporary artifacts, and scratch work (gitignored).

## Python and Tooling

- Target Python version: 3.11.
- Use `uv` for dependency management and command execution.
- Use `ruff` for linting and formatting.
- Use `pytest` for tests.

Typical commands:

- `uv sync`
- `uv run ruff check .`
- `uv run ruff format .`
- `uv run pytest`

Test commands:

- Run all tests: `uv run pytest`
- Run one test file: `uv run pytest tests/test_arithmetic.py`
- Run one test by name: `uv run pytest -k test_div`

## C++ Guidelines and Tooling

- Prefer modern C++ style (C++17 or newer).
- Prefer RAII and standard library ownership types (`std::unique_ptr`, `std::shared_ptr`) over raw `new`/`delete`.
- Prefer `std::vector`, `std::array`, and standard algorithms over manual memory management.
- Use `const` correctness, references where appropriate, and explicit numeric types for reproducibility.
- Keep numerics and model logic in reusable library code; keep executables as thin orchestration layers.
- Use repository `.clang-format` and `.clangd` settings for consistent formatting and editor tooling.

CMake/CTest workflow:

- Configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
- Build: `cmake --build build -j`
- Run tests: `ctest --test-dir build --output-on-failure`
- Run one test: `ctest --test-dir build -R <test_name> --output-on-failure`

## Notebook Policy (Marimo)

- Prefer marimo for all new interactive research workflows.
- Place marimo notebooks in `notebooks/`.
- Do not add new Jupyter notebooks (`.ipynb`) unless there is a compelling external requirement.
- Keep notebooks thin: import reusable logic from `src/sample/`.

## Validation and Reproducibility

- For physics-impacting changes, include at least one validation signal:
  - automated test,
  - numerical invariant check,
  - or documented benchmark comparison.
- Record seeds, key parameters, and command lines when producing publishable outputs.
- Avoid silent changes to assumptions, constants, units, or model behavior.

## Branch and Pages Policy

- `main`: stable branch.
- `develop`: active integration branch.
- `feature/*`: short-lived feature work.
- Prefer not to push directly to `main`.
- Do most active development in `develop`; merge to `main` when changes are stable.
- Prefer a normal merge commit when merging `develop` into `main`.
- Merge feature branches into `develop` first, then flow those changes to `main` via `develop`.
- For `feature/*` into `develop`, squash-and-merge is preferred to avoid noisy micro-commits.

GitHub Pages publishing policy:

- Publish only `main` and `develop`.
- Use paths:
  - `/main/` from `main`
  - `/develop/` from `develop`
- Do not auto-publish `feature/*` branches.
- Notebook export mode defaults to static HTML.
- List notebooks in `notebooks/publish.toml` under `wasm` to publish those as interactive WASM pages.

## Data and Artifact Handling

- Keep large files, generated figures, caches, and temporary data in `work/`.
- Commit only small fixtures needed for tests or reproducible examples.
- Prefer deterministic scripts over manual steps when generating artifacts.

## Agent Behavior Expectations

- Explain what changed and why, especially for scientific logic.
- Prefer conservative edits that preserve existing behavior unless a change is intentional.
- When changing scientific behavior, include validation evidence in the change description.
