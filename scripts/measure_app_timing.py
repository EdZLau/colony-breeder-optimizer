"""Measure synthetic Streamlit AppTest ranking and display-only rerun timings.

All colony rows are synthetic. Results are single-machine measurements, not
performance guarantees.
"""

import csv
from io import StringIO
from pathlib import Path
from time import perf_counter

from streamlit.testing.v1 import AppTest

APP_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "colony_breeder"
    / "app"
    / "streamlit_app.py"
)
RUN_TIMEOUT_SECONDS = 180
REPETITIONS = 3
ANIMALS_PER_SEX = 200


def _csv_bytes(headers, rows):
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\r\n")
    writer.writerow(headers)
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def _scenario_a_csv():
    headers = ["animal_id", "sex", "cage", "mother_id", "father_id", "A"]
    rows = []
    for sex, prefix in (("M", "M"), ("F", "F")):
        for index in range(ANIMALS_PER_SEX):
            rows.append(
                (f"{prefix}{index:03}", sex, "", "", "", "A/a")
            )
    return _csv_bytes(headers, rows)


def _scenario_b_csv():
    headers = [
        "animal_id",
        "sex",
        "cage",
        "mother_id",
        "father_id",
        "U1",
        "U2",
        "L1",
        "L2",
        "phase_L1_L2",
    ]
    rows = []
    for sex, prefix in (("M", "M"), ("F", "F")):
        for index in range(ANIMALS_PER_SEX):
            missing_target = sex == "M" and index < 20
            rows.append(
                (
                    f"{prefix}{index:03}",
                    sex,
                    "",
                    "",
                    "",
                    "" if missing_target else "A/a",
                    "B/b",
                    "D/d",
                    "E/e",
                    "D|E/d|e",
                )
            )
    return _csv_bytes(headers, rows)


def _widget(app, widgets, label):
    return next(widget for widget in widgets if widget.label == label)


def _select_pool_genes(app, genes):
    for gene in genes:
        _widget(
            app,
            app.multiselect,
            "Gene columns to use in settings",
        ).select(gene).run(timeout=RUN_TIMEOUT_SECONDS)


def _select_target_gene(app, gene):
    _widget(app, app.multiselect, "Target genes").select(gene).run(
        timeout=RUN_TIMEOUT_SECONDS
    )


def _set_target_text(app, gene, genotype):
    widget = next(
        item
        for item in app.text_area
        if item.label.startswith(f"Target genotypes for {gene}")
    )
    widget.set_value(genotype).run(timeout=RUN_TIMEOUT_SECONDS)


def _prepare_app(scenario):
    if scenario == "A":
        content = _scenario_a_csv()
        pool_genes = ("A",)
        targets = (("A", "A/a"),)
    else:
        content = _scenario_b_csv()
        pool_genes = ("U1", "U2", "L1", "L2")
        targets = (("U1", "A/a"), ("U2", "B/b"), ("L1", "D/d"))

    app = AppTest.from_file(str(APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        f"synthetic-scenario-{scenario}.csv",
        content,
        "text/csv",
    ).run(timeout=RUN_TIMEOUT_SECONDS)
    _widget(app, app.button, "Load colony").click().run(
        timeout=RUN_TIMEOUT_SECONDS
    )
    _select_pool_genes(app, pool_genes)
    _widget(app, app.button, "Load colony").click().run(
        timeout=RUN_TIMEOUT_SECONDS
    )

    for gene, genotype in targets:
        _select_target_gene(app, gene)
        _set_target_text(app, gene, genotype)

    if scenario == "B":
        _widget(app, app.number_input, "Linked gene pairs (up to 3)").set_value(
            1
        ).run(timeout=RUN_TIMEOUT_SECONDS)
        _widget(app, app.selectbox, "Gene 1").select("L1").run(
            timeout=RUN_TIMEOUT_SECONDS
        )
        _widget(app, app.selectbox, "Gene 2").select("L2").run(
            timeout=RUN_TIMEOUT_SECONDS
        )
        _widget(app, app.text_input, "Recombination fraction r").set_value(
            "1/10"
        ).run(timeout=RUN_TIMEOUT_SECONDS)
        _widget(app, app.selectbox, "Phase column").select("phase_L1_L2").run(
            timeout=RUN_TIMEOUT_SECONDS
        )

    if app.exception:
        raise RuntimeError(f"Scenario {scenario} setup raised: {app.exception}")
    if not any(
        element.value == "Settings valid."
        for element in app.success
    ):
        errors = [element.value for element in app.error]
        raise RuntimeError(
            f"Scenario {scenario} settings are invalid: {errors}"
        )
    return app


def _time_run_and_display_rerun(scenario):
    app = _prepare_app(scenario)
    started = perf_counter()
    _widget(app, app.button, "Run ranking").click().run(
        timeout=RUN_TIMEOUT_SECONDS
    )
    ranking_seconds = perf_counter() - started

    result = app.session_state["run_result"]
    if result.error_message is not None:
        raise RuntimeError(
            f"Scenario {scenario} ranking failed: {result.error_message}"
        )
    if "Ranking results" not in [item.value for item in app.subheader]:
        raise RuntimeError(f"Scenario {scenario} did not render results.")

    started = perf_counter()
    _widget(app, app.number_input, "Top N pairs to show").set_value(5).run(
        timeout=RUN_TIMEOUT_SECONDS
    )
    top_rerun_seconds = perf_counter() - started
    if "Ranking results" not in [item.value for item in app.subheader]:
        raise RuntimeError(
            f"Scenario {scenario} top-N-only rerun hid the result."
        )

    return len(result.ranked), len(result.skipped), ranking_seconds, top_rerun_seconds


def main():
    print("Synthetic AppTest timing; tracemalloc is not enabled.")
    print("Times are single-machine measurements in seconds.")
    for scenario, description in (
        ("A", "one gene, uniform A/a, target A/a, no skips"),
        (
            "B",
            "two unlinked genes + one linked pair; U1 blank for 20 males "
            "(4,000 missing-genotype skips)",
        ),
    ):
        print(f"\nScenario {scenario}: {description}")
        for repetition in range(1, REPETITIONS + 1):
            ranked, skipped, run_seconds, rerun_seconds = (
                _time_run_and_display_rerun(scenario)
            )
            print(
                f"  run {repetition}: ranked={ranked}, skipped={skipped}, "
                f"Run ranking click={run_seconds:.3f}s, "
                f"top-N-only rerun={rerun_seconds:.3f}s"
            )


if __name__ == "__main__":
    main()
