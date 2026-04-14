import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import sample

    return mo, sample


@app.cell
def _(mo):
    from textwrap import dedent

    mo.md(
        dedent(
            """
            # Marimo Notebook Demo (Static)

            This notebook is intended for local runs and static publication.
            It can import the local `sample` module from `src/sample/` and run
            simple checks for project utility functions.
            """
        )
    )
    return


@app.cell
def _(sample):
    arithmetic_results = {
        "add(2, 3)": sample.add(2, 3),
        "sub(7, 4)": sample.sub(7, 4),
        "mul(6, 5)": sample.mul(6, 5),
        "div(8, 2)": sample.div(8, 2),
    }

    assert arithmetic_results["add(2, 3)"] == 5
    assert arithmetic_results["sub(7, 4)"] == 3
    assert arithmetic_results["mul(6, 5)"] == 30
    assert arithmetic_results["div(8, 2)"] == 4

    div_zero_ok = False
    try:
        sample.div(1, 0)
    except ValueError:
        div_zero_ok = True
    assert div_zero_ok

    return (arithmetic_results,)


@app.cell
def _(arithmetic_results, mo):
    rows = "\n".join(f"- `{expr}` -> `{value}`" for expr, value in arithmetic_results.items())
    mo.md(
        "## sample module checks\n\n"
        "This notebook imports the local `sample` module and verifies arithmetic behavior.\n\n"
        f"{rows}\n\n"
        "- `div(1, 0)` raises `ValueError` as expected"
    )
    return


if __name__ == "__main__":
    app.run()
