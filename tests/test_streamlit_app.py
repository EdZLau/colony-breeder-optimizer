from pathlib import Path

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
_APP_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "colony_breeder"
    / "app"
    / "streamlit_app.py"
)


def _rendered_markdown(app):
    return [element.value for element in app.markdown]


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
