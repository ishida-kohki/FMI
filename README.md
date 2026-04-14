# Skelton Project

This repository is a starter template for a scientific research project.

## What you can do here

- Build reusable code for models and analysis
- Run interactive notebooks with marimo
- Track tests and documentation as the project grows
- Publish notebook results to GitHub Pages (`/main/` and `/develop/`)

## Published pages

- Landing page: <https://amanotk.github.io/skelton-project/>
- Stable notebooks (`main`): <https://amanotk.github.io/skelton-project/main/>
- Development notebooks (`develop`): <https://amanotk.github.io/skelton-project/develop/>

## Directory structure

- `src/`: reusable Python code
- `notebooks/`: interactive marimo notebooks
- `tests/`: automated tests
- `docs/`: project documentation
- `work/`: large files and scratch outputs (not committed)

## Quick start

1. Install dependencies:

```bash
uv sync
```

2. Run tests:

```bash
uv run pytest
```

3. Open a notebook:

```bash
uv run marimo run notebooks/demo_static.py
```

WASM-focused example notebook:

```bash
uv run marimo run notebooks/demo_wasm.py
```

## Marimo notebook workflow

- Run in browser: `uv run marimo run notebooks/<notebook-name>.py`
- Edit notebook: `uv run marimo edit notebooks/<notebook-name>.py`
- Script-mode check: `uv run notebooks/<notebook-name>.py`
- Notebook lint check: `uvx marimo check notebooks/<notebook-name>.py`

## GitHub Pages export mode

- Default behavior: notebooks are published as static HTML.
- To publish a notebook as interactive WASM, list it in `notebooks/publish.toml`.

Example:

```toml
wasm = [
  "notebooks/demo_wasm.py",
]
```

`demo_static.py` imports the local `sample` module; `demo_wasm.py` avoids local imports so it works in browser WASM mode.

## Recommended Git branching strategy

- `main`: stable work (notebooks published to `/main/` on GitHub Pages)
- `develop`: active development (notebooks published to `/develop/` on GitHub Pages)
- `feature/*`: short-lived feature branches

### Development on `develop`
- Commit primarily to `develop` for ongoing work
- Merge `develop` into `main` when stable and ready for release (usually with a normal merge commit)

### Feature branches
- Create `feature/*` branches for specific features or experiments
- Merge `feature/*` branches back into `develop` when ready (usually with squash merge)
