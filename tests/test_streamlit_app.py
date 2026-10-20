from pathlib import Path
import re

import pytest

pytest.importorskip("streamlit")

from streamlit.testing.v1 import AppTest

_CSV_CONTENT = (
    b"Animal ID,Sex,Cage,Dam,Sire,Cre\n"
    b"M01,M,1,,,Cre/+\n"
    b"F01,,1,,,Cre/+\n"
)
_CSV_WITHOUT_SEX = (
    b"Animal ID,Cage,Dam,Sire,Cre\n"
    b"M01,1,,,Cre/+\n"
)
_SETTINGS_CSV = (
    b"Animal ID,Sex,Cage,Dam,Sire,Cre,Xgene,A,B,AB_phase\n"
    b"M01,M,C1,D1,S1,Cre/+,A/Y,A/a,B/b,A|B/a|b\n"
    b"F01,,C2,D1,S2,Cre/+,A/a,A/a,B/b,A|B/a|b\n"
)
_APP_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "colony_breeder"
    / "app"
    / "streamlit_app.py"
)


def _rendered_markdown(app):
    return [element.value for element in app.markdown]


def _loaded_settings_app(gene_columns):
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "settings.csv",
        _SETTINGS_CSV,
        "text/csv",
    ).run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )
    pool = next(
        widget
        for widget in app.multiselect
        if widget.label == "Gene columns to inspect (not parsed in this step)"
    )
    for gene in gene_columns:
        pool.select(gene).run(timeout=30)
        pool = next(
            widget
            for widget in app.multiselect
            if widget.label == "Gene columns to inspect (not parsed in this step)"
        )
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )
    assert not app.exception
    return app


def _set_valid_target(app, gene, genotype):
    target_genes = next(
        widget for widget in app.multiselect if widget.label == "Target genes"
    )
    target_genes.select(gene).run(timeout=30)
    target_text = next(
        widget
        for widget in app.text_area
        if widget.label.startswith(f"Target genotypes for {gene}")
    )
    target_text.set_value(genotype).run(timeout=30)
    assert not app.exception


def test_streamlit_page_renders_upload_widget():
    app = AppTest.from_file(str(_APP_PATH)).run()

    assert not app.exception
    assert len(app.get("file_uploader")) == 1


def test_upload_and_load_shows_animal_and_unknown_sex_counts():
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "tiny.csv",
        _CSV_CONTENT,
        "text/csv",
    ).run(timeout=30)
    assert not app.exception

    load_button = next(button for button in app.button if button.label == "Load colony")
    assert not load_button.disabled
    load_button.click().run(timeout=30)

    assert not app.exception
    rendered = " ".join(_rendered_markdown(app))
    assert "Loaded animals: 2" in rendered
    assert "Animals with unknown sex: 1" in rendered
    assert sum(
        "Genotypes have NOT been loaded or validated" in element.value
        for element in app.caption
    ) == 1
    assert sum(
        "Genotypes have NOT been loaded or validated" in element.value
        for element in app.info
    ) == 1


def test_changing_header_row_clears_loaded_result():
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "tiny.csv",
        _CSV_CONTENT,
        "text/csv",
    ).run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )

    assert "Load result" in [element.value for element in app.subheader]
    app.number_input[0].set_value(2).run(timeout=30)

    assert not app.exception
    assert "Load result" not in [element.value for element in app.subheader]


def test_load_is_disabled_when_no_sex_column_is_available():
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "no-sex.csv",
        _CSV_WITHOUT_SEX,
        "text/csv",
    ).run(timeout=30)

    assert not app.exception
    load_button = next(button for button in app.button if button.label == "Load colony")
    assert load_button.disabled
    assert any(
        "Load disabled: choose a sex column" in element.value
        for element in app.caption
    )


def test_load_is_disabled_when_sex_column_has_too_many_distinct_values():
    csv_content = (
        b"Animal ID,Sex,Cage,Dam,Sire,Cre\n"
        + b"".join(
            f"M{index:02},value{index:02},1,,,Cre/+\n".encode()
            for index in range(21)
        )
    )
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "many-sex-values.csv",
        csv_content,
        "text/csv",
    ).run(timeout=30)

    assert not app.exception
    load_button = next(button for button in app.button if button.label == "Load colony")
    assert load_button.disabled
    assert any(
        "Load disabled:" in element.value
        and "more than the limit" in element.value
        for element in app.caption
    )


def test_settings_controls_validate_live_without_running_ranking():
    app = AppTest.from_file(str(_APP_PATH)).run()
    app.get("file_uploader")[0].upload(
        "tiny.csv",
        _CSV_CONTENT,
        "text/csv",
    ).run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )

    gene_pool = next(
        widget
        for widget in app.multiselect
        if widget.label == "Gene columns to inspect (not parsed in this step)"
    )
    gene_pool.select("Cre").run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )

    assert not app.exception
    target_genes = next(
        widget for widget in app.multiselect if widget.label == "Target genes"
    )
    target_genes.select("Cre").run(timeout=30)

    assert any("Target genotype list" in element.value for element in app.markdown)
    target_text = next(
        widget for widget in app.text_area if widget.label.startswith("Target genotypes for Cre")
    )
    target_text.set_value("Cre/+").run(timeout=30)

    assert not app.exception
    assert any(element.value == "Settings valid." for element in app.success)
    assert any("ranking is not available yet" in element.value for element in app.caption)
    assert all(button.label != "Run ranking" for button in app.button)
    assert "Ranked results" not in [element.value for element in app.subheader]


def test_settings_signature_uses_upload_inputs_but_ignores_top():
    app = _loaded_settings_app(["Cre", "Xgene"])
    _set_valid_target(app, "Cre", "Cre/+")
    baseline = app.session_state["settings_signature"]

    top = next(widget for widget in app.number_input if widget.label == "Top N pairs to show")
    top.set_value(25).run(timeout=30)
    assert not app.exception
    assert app.session_state["settings_signature"] == baseline

    app = _loaded_settings_app(["Cre", "Xgene"])
    _set_valid_target(app, "Cre", "Cre/+")
    baseline = app.session_state["settings_signature"]
    sex_widget = next(widget for widget in app.selectbox if widget.label == "(blank)")
    sex_widget.select("F").run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )
    assert not app.exception
    assert app.session_state["settings_signature"] != baseline

    app = _loaded_settings_app(["Cre", "Xgene"])
    _set_valid_target(app, "Cre", "Cre/+")
    baseline = app.session_state["settings_signature"]
    cage_mapping = next(widget for widget in app.selectbox if widget.label == "Cage")
    cage_mapping.select("(none)").run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )
    assert not app.exception
    assert app.session_state["settings_signature"] != baseline

    app = _loaded_settings_app(["Cre", "Xgene"])
    _set_valid_target(app, "Cre", "Cre/+")
    baseline = app.session_state["settings_signature"]
    pool = next(
        widget
        for widget in app.multiselect
        if widget.label == "Gene columns to inspect (not parsed in this step)"
    )
    pool.select("A").run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )
    assert not app.exception
    assert app.session_state["settings_signature"] != baseline


def test_invalid_genotype_shows_bulleted_errors_without_valid_state():
    app = _loaded_settings_app(["Cre"])
    _set_valid_target(app, "Cre", "Cre/")

    assert not app.exception
    assert any(
        "- Invalid target genotype" in element.value
        for element in app.markdown
    )
    assert not any(element.value == "Settings valid." for element in app.success)


@pytest.mark.parametrize(
    ("genotype", "valid"),
    [
        ("A", False),
        ("A/B/C", False),
        ("/Y", False),
        ("A/A", True),
        ("A/Y", True),
    ],
)
def test_x_linked_target_genotypes_follow_core_parser(genotype, valid):
    app = _loaded_settings_app(["Xgene"])
    x_linked = next(widget for widget in app.multiselect if widget.label == "X-linked genes")
    x_linked.select("Xgene").run(timeout=30)
    _set_valid_target(app, "Xgene", genotype)

    assert not app.exception
    assert any(element.value == "Settings valid." for element in app.success) is valid
    assert any("- Invalid target genotype" in element.value for element in app.markdown) is not valid


@pytest.mark.parametrize(
    ("recombination", "valid"),
    [
        ("0.6", False),
        ("abc", False),
        ("1/0", False),
        ("1/10", True),
    ],
)
def test_linked_pair_recombination_text_is_validated(recombination, valid):
    app = _loaded_settings_app(["A", "B"])
    _set_valid_target(app, "A", "A/a")
    count = next(
        widget for widget in app.number_input
        if widget.label == "Linked gene pairs (up to 3)"
    )
    count.set_value(1).run(timeout=30)

    next(widget for widget in app.selectbox if widget.label == "Gene 1").select("A").run(
        timeout=30
    )
    next(widget for widget in app.selectbox if widget.label == "Gene 2").select("B").run(
        timeout=30
    )
    next(widget for widget in app.selectbox if widget.label == "Phase column").select(
        "AB_phase"
    ).run(timeout=30)
    next(widget for widget in app.text_input if widget.label == "Recombination fraction r").set_value(
        recombination
    ).run(timeout=30)

    assert not app.exception
    assert any(element.value == "Settings valid." for element in app.success) is valid
    assert any("- Invalid recombination fraction" in element.value for element in app.markdown) is not valid


def test_three_blank_linked_pair_rows_do_not_add_validation_errors():
    app = _loaded_settings_app(["Cre"])
    _set_valid_target(app, "Cre", "Cre/+")
    count = next(
        widget for widget in app.number_input
        if widget.label == "Linked gene pairs (up to 3)"
    )
    count.set_value(3).run(timeout=30)

    assert not app.exception
    assert sum(widget.label == "Gene 1" for widget in app.selectbox) == 3
    assert any(element.value == "Settings valid." for element in app.success)


def test_linked_pair_gene_options_use_the_eligible_gene_pool():
    app = _loaded_settings_app(["A", "B", "Cre"])
    _set_valid_target(app, "Cre", "Cre/+")
    count = next(
        widget for widget in app.number_input
        if widget.label == "Linked gene pairs (up to 3)"
    )
    count.set_value(1).run(timeout=30)

    gene1 = next(widget for widget in app.selectbox if widget.label == "Gene 1")
    assert gene1.options == ["(none)", "A", "B", "Cre"]


def test_removing_target_gene_removes_text_area_and_setting():
    app = _loaded_settings_app(["Cre"])
    _set_valid_target(app, "Cre", "Cre/+")

    target_genes = next(widget for widget in app.multiselect if widget.label == "Target genes")
    target_genes.unselect("Cre").run(timeout=30)

    assert not app.exception
    assert not any(
        widget.label.startswith("Target genotypes for Cre")
        for widget in app.text_area
    )
    assert dict(app.session_state["run_settings_draft"].targets) == {}


def test_removing_gene_from_pool_cleans_stale_linked_selectbox_state():
    app = _loaded_settings_app(["A", "B", "Cre"])
    _set_valid_target(app, "Cre", "Cre/+")
    count = next(
        widget for widget in app.number_input
        if widget.label == "Linked gene pairs (up to 3)"
    )
    count.set_value(1).run(timeout=30)
    next(widget for widget in app.selectbox if widget.label == "Gene 1").select("A").run(
        timeout=30
    )
    next(widget for widget in app.selectbox if widget.label == "Gene 2").select("B").run(
        timeout=30
    )

    pool = next(
        widget
        for widget in app.multiselect
        if widget.label == "Gene columns to inspect (not parsed in this step)"
    )
    pool.unselect("A").run(timeout=30)
    next(button for button in app.button if button.label == "Load colony").click().run(
        timeout=30
    )

    assert not app.exception
    assert next(
        widget for widget in app.number_input
        if widget.label == "Linked gene pairs (up to 3)"
    ).value == 0
    assert app.session_state["run_settings_draft"].linked_pairs == ()
    assert "A" not in app.session_state["run_settings_draft"].targets


def test_page_source_does_not_import_or_call_core_ranking_or_colony_loader():
    source = _APP_PATH.read_text(encoding="utf-8")

    assert "run_ranking" not in source
    assert "rank_pairs" not in source
    assert not re.search(r"\bload_colony(?:_rows)?\s*\(", source)
