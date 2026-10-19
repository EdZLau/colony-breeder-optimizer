from fractions import Fraction

import pytest

from colony_breeder.app.run_settings import (
    build_reload_args,
    build_run_settings,
    parse_recombination_fraction,
    validate_run_settings,
)
from colony_breeder.app.ui_helpers import build_loader_args
from colony_breeder.core.colony import load_colony_rows


FIELDNAMES = ["animal_id", "sex", "cage", "mother_id", "father_id", "A", "B", "AB_phase"]
MAPPING = {
    "animal_id": "animal_id",
    "sex": "sex",
    "cage": "cage",
    "mother_id": "mother_id",
    "father_id": "father_id",
}
GENE_POOL = ["A", "B"]


def _settings(**overrides):
    values = {
        "targets": {"A": ["A/a"]},
        "lethals": {},
        "x_linked_genes": (),
        "linked_pairs": (),
        "exclude_related": False,
        "top": None,
    }
    values.update(overrides)
    return build_run_settings(**values)


def _errors(settings, gene_pool=GENE_POOL, fieldnames=FIELDNAMES, mapping=MAPPING):
    return validate_run_settings(settings, gene_pool, fieldnames, mapping)


def test_build_run_settings_strips_and_stably_deduplicates_values():
    settings = build_run_settings(
        targets={" A ": [" A/a ", "A/a", "a/a", " a/a "]},
        lethals={" B ": " B/B "},
        x_linked_genes=[" X ", "X", "Y"],
        linked_pairs=[(" A ", " B ", " 1/10 ", " A_B_phase ")],
        exclude_related=True,
        top=5,
    )

    assert dict(settings.targets) == {"A": ("A/a", "a/a")}
    assert dict(settings.lethals) == {"B": ("B/B",)}
    assert settings.x_linked_genes == ("X", "Y")
    assert settings.linked_pairs == (("A", "B", "1/10", "A_B_phase"),)
    assert settings.exclude_related is True
    assert settings.top == 5
    with pytest.raises(TypeError):
        settings.targets["A"] = ("A/A",)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("1/10", Fraction(1, 10)),
        ("0.1", Fraction(1, 10)),
        ("1/2", Fraction(1, 2)),
        ("0", Fraction(0)),
    ],
)
def test_recombination_text_is_parsed_as_exact_fraction(text, expected):
    assert parse_recombination_fraction(text) == expected
    assert _errors(
        _settings(linked_pairs=[("A", "B", text, "AB_phase")])
    ) == []


@pytest.mark.parametrize("text", ["0.6", "abc"])
def test_out_of_range_or_malformed_recombination_text_is_a_validation_error(text):
    errors = _errors(
        _settings(linked_pairs=[("A", "B", text, "AB_phase")])
    )
    assert len(errors) == 1
    assert "recombination fraction" in errors[0]


@pytest.mark.parametrize("r_value", [True, 0.1])
def test_non_text_recombination_input_is_rejected(r_value):
    with pytest.raises(ValueError, match="must be text"):
        build_run_settings(
            targets={"A": "A/a"},
            linked_pairs=[("A", "B", r_value, "AB_phase")],
        )


def test_settings_require_at_least_one_target():
    assert _errors(_settings(targets={})) == [
        "Choose at least one target gene and genotype."
    ]


def test_target_gene_must_be_in_selected_gene_pool():
    errors = _errors(_settings(targets={"C": "C/c"}))
    assert len(errors) == 1
    assert "Target gene 'C'" in errors[0]


def test_lethal_gene_must_be_in_selected_gene_pool():
    errors = _errors(_settings(lethals={"C": "C/C"}))
    assert len(errors) == 1
    assert "Lethal gene 'C'" in errors[0]


def test_linked_genes_must_be_in_selected_gene_pool():
    errors = _errors(
        _settings(linked_pairs=[("A", "C", "1/10", "AB_phase")])
    )
    assert errors == [
        "Linked gene 'C' is not in the selected gene columns."
    ]


def test_x_linked_gene_must_be_in_selected_gene_pool():
    errors = _errors(_settings(x_linked_genes=["X"]))
    assert errors == ["X-linked gene 'X' is not in the selected gene columns."]


def test_autosomal_target_genotype_is_validated_by_core_parser():
    errors = _errors(_settings(targets={"A": "A++"}))
    assert len(errors) == 1
    assert "Invalid target genotype 'A++'" in errors[0]


def test_autosomal_lethal_genotype_is_validated_by_core_parser():
    errors = _errors(_settings(lethals={"A": "A++"}))
    assert len(errors) == 1
    assert "Invalid lethal genotype 'A++'" in errors[0]


def test_x_linked_genotype_requires_slash_form_and_male_y():
    errors = _errors(
        _settings(targets={"A": "A/Y"}, x_linked_genes=["A"], lethals={"A": "A/y"})
    )
    assert errors == []

    malformed = _errors(
        _settings(targets={"A": "A/Y"}, x_linked_genes=["A"], lethals={"A": "A/Y/Y"})
    )
    assert len(malformed) == 1
    assert "Invalid lethal genotype" in malformed[0]


def test_target_genotype_list_cannot_be_empty():
    errors = _errors(_settings(targets={"A": []}))
    assert errors == ["Target genotype list for gene 'A' cannot be empty."]


def test_linked_pair_cannot_repeat_gene_within_pair():
    errors = _errors(
        _settings(linked_pairs=[("A", "A", "1/10", "AB_phase")])
    )
    assert any("same gene twice" in error for error in errors)


def test_linked_gene_cannot_repeat_across_pairs():
    errors = _errors(
        _settings(
            linked_pairs=[
                ("A", "B", "1/10", "AB_phase"),
                ("B", "A", "1/10", "BA_phase"),
            ]
        )
    )
    assert sum("repeated within or across" in error for error in errors) == 2


def test_linked_gene_cannot_also_be_x_linked():
    errors = _errors(
        _settings(
            x_linked_genes=["A"],
            linked_pairs=[("A", "B", "1/10", "AB_phase")],
        )
    )
    assert any("cannot include X-linked gene 'A'" in error for error in errors)


def test_linked_pair_requires_existing_phase_column():
    errors = _errors(
        _settings(linked_pairs=[("A", "B", "1/10", "Missing_phase")])
    )
    assert len(errors) == 1
    assert "not present in the uploaded headers" in errors[0]


def test_phase_column_cannot_be_mapped_to_colony_field():
    mapping = dict(MAPPING, cage="AB_phase")
    errors = _errors(
        _settings(linked_pairs=[("A", "B", "1/10", "AB_phase")]),
        mapping=mapping,
    )
    assert any("cannot be mapped to a colony field" in error for error in errors)


def test_phase_column_cannot_be_a_gene_column():
    errors = _errors(
        _settings(linked_pairs=[("A", "B", "1/10", "A")])
    )
    assert any("cannot also be a selected gene column" in error for error in errors)


def test_phase_column_cannot_be_shared_between_linked_pairs():
    errors = _errors(
        _settings(
            linked_pairs=[
                ("A", "B", "1/10", "AB_phase"),
                ("B", "A", "1/10", "AB_phase"),
            ]
        )
    )
    assert any("cannot be shared between linked pairs" in error for error in errors)


@pytest.mark.parametrize("top", [0, -1, True, 1.5])
def test_top_must_be_a_positive_integer_when_set(top):
    errors = _errors(_settings(top=top))
    assert errors == ["Top N must be an integer greater than or equal to 1."]


def test_reload_arguments_use_ordered_gene_union_and_exact_recombination():
    settings = _settings(
        targets={"B": "B/b", "A": "A/a"},
        lethals={"B": "B/B", "C": "C/C"},
        x_linked_genes=["D"],
        linked_pairs=[("E", "A", "0.1", "EA_phase")],
    )

    args = build_reload_args({"first_data_row": 2}, settings)

    assert args["gene_columns"] == ["B", "A", "C", "D", "E"]
    assert args["x_linked_genes"] == ["D"]
    assert args["linked_pairs"] == [
        ("E", "A", Fraction(1, 10), "EA_phase")
    ]
    assert "EA_phase" not in args["gene_columns"]
    assert args["allow_unknown_sex"] is True
    assert args["first_data_row"] == 2


def test_reload_arguments_load_a_linked_pair_from_tiny_fake_rows():
    rows = [
        {
            "animal_id": "M01",
            "sex": "M",
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "A": "A/a",
            "B": "B/b",
            "AB_phase": "A|B/a|b",
        }
    ]
    upload_args = build_loader_args(FIELDNAMES, rows, MAPPING, first_data_row=2)
    settings = _settings(linked_pairs=[("A", "B", "1/10", "AB_phase")])

    animals = load_colony_rows(
        **build_reload_args(upload_args, settings)
    )

    assert animals[0].genotypes == {"A": "A/a", "B": "B/b"}
    assert animals[0].phases == {("A", "B"): "A|B/a|b"}
