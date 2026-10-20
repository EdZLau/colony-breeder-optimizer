from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction

import pytest

from colony_breeder.app.run_ranking import (
    ExcludedAnimalView,
    RunResult,
    SkippedPairView,
    result_is_current,
    run_ranking,
    skipped_counts,
    visible_ranked,
)
from colony_breeder.app.run_settings import build_run_settings, run_signature
from colony_breeder.app.ui_helpers import build_loader_args


MAPPING = {
    "animal_id": "animal_id",
    "sex": "sex",
    "cage": "cage",
    "mother_id": "mother_id",
    "father_id": "father_id",
}
BASE_FIELDS = ["animal_id", "sex", "cage", "mother_id", "father_id"]


def _row(animal_id, sex, **values):
    row = {
        "animal_id": animal_id,
        "sex": sex,
        "cage": f"C{animal_id}",
        "mother_id": "",
        "father_id": "",
    }
    row.update(values)
    return row


def _run(rows, settings, genes, extra_fields=()):
    fieldnames = [*BASE_FIELDS, *genes, *extra_fields]
    upload_args = build_loader_args(
        fieldnames,
        rows,
        MAPPING,
        first_data_row=2,
    )
    return run_ranking(
        upload_args,
        settings,
        gene_pool=genes,
        fieldnames=fieldnames,
        column_mapping=MAPPING,
        base_signature=("upload", "tiny"),
    )


def test_autosomal_a_a_cross_targeting_recessive_homozygote_is_one_quarter():
    rows = [
        _row("M01", "M", A="A/a"),
        _row("F01", "F", A="A/a"),
    ]
    result = _run(rows, build_run_settings(targets={"A": "a/a"}), ["A"])

    assert result.error_message is None
    assert result.ranked[0].probability == Fraction(1, 4)
    assert result.ranked[0].per_gene == (("A", Fraction(1, 4)),)


def test_invalid_settings_return_error_message_without_ranked_pairs():
    rows = [_row("M01", "M", A="A/a"), _row("F01", "F", A="A/a")]
    result = _run(rows, build_run_settings(targets={"A": "A++"}), ["A"])

    assert result.error_message is not None
    assert "Invalid target genotype" in result.error_message
    assert result.ranked == ()
    assert result.skipped == ()


def test_malformed_phase_is_returned_as_error_not_raised():
    rows = [
        _row("M01", "M", A="A/a", B="B/b", AB_phase="bad phase"),
        _row("F01", "F", A="A/a", B="B/b", AB_phase="A|B/a|b"),
    ]
    settings = build_run_settings(
        targets={"A": "A/a"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
    )

    result = _run(rows, settings, ["A", "B"], ["AB_phase"])

    assert result.error_message is not None
    assert "Invalid phase" in result.error_message
    assert result.ranked == ()


def test_blank_target_genotype_is_reported_as_missing_genotype_skip():
    rows = [
        _row("M01", "M", A="A/a"),
        _row("F01", "F", A=""),
    ]
    result = _run(rows, build_run_settings(targets={"A": "A/a"}), ["A"])

    assert result.ranked == ()
    assert len(result.skipped) == 1
    assert result.skipped[0].kind == "missing_genotype"
    assert "female F01" in result.skipped[0].reason
    assert skipped_counts(result)["missing_genotype"] == 1


def test_unknown_sex_is_reported_as_excluded_and_not_paired():
    rows = [
        _row("M01", "M", A="A/a"),
        _row("F01", "F", A="A/a"),
        _row("U01", "", A="A/a"),
    ]
    result = _run(rows, build_run_settings(targets={"A": "A/a"}), ["A"])

    assert len(result.ranked) == 1
    assert {(animal.animal_id, animal.reason) for animal in result.excluded_animals} == {
        (
            "U01",
            "sex is blank or unrecognized, so it can't be assigned a male or female role",
        )
    }
    assert result.notes == (
        "Some animals were excluded because sex is blank or unrecognized, "
        "so they can't be assigned a male or female role.",
    )


def test_pedigree_warning_text_is_preserved_from_core():
    rows = [
        _row("M01", "M", A="A/a", mother_id="RETIRED_DAM"),
        _row("F01", "F", A="A/a"),
    ]
    result = _run(rows, build_run_settings(targets={"A": "A/a"}), ["A"])

    assert result.pedigree_warnings == (
        (
            "info",
            "Animal 'M01' has mother_id 'RETIRED_DAM', which is not in the colony.",
        ),
    )


def test_x_linked_run_uses_core_probability():
    rows = [
        _row("M01", "M", X="A/Y"),
        _row("F01", "F", X="A/a"),
    ]
    settings = build_run_settings(
        targets={"X": "A/Y"},
        x_linked_genes=["X"],
    )
    result = _run(rows, settings, ["X"])

    assert result.error_message is None
    assert result.ranked[0].probability == Fraction(1, 4)
    assert result.ranked[0].per_gene == (("X", Fraction(1, 4)),)


def test_linked_run_with_one_tenth_recombination_is_exact():
    rows = [
        _row("M01", "M", A="A/a", B="B/b", AB_phase="A|B/a|b"),
        _row("F01", "F", A="A/a", B="B/b", AB_phase="A|B/a|b"),
    ]
    settings = build_run_settings(
        targets={"A": "A/A", "B": "B/B"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
    )
    result = _run(rows, settings, ["A", "B"], ["AB_phase"])

    assert result.error_message is None
    assert result.ranked[0].probability == Fraction(81, 400)
    assert result.ranked[0].linked_pair_probabilities == (
        (("A", "B"), Fraction(81, 400)),
    )


def test_exclude_related_keeps_core_flags_on_skipped_pair():
    rows = [
        _row("M01", "M", A="A/a", mother_id="MOTHER", father_id="FATHER"),
        _row("F01", "F", A="A/a", mother_id="MOTHER", father_id="FATHER"),
    ]
    result = _run(
        rows,
        build_run_settings(targets={"A": "A/a"}, exclude_related=True),
        ["A"],
    )

    assert result.ranked == ()
    assert len(result.skipped) == 1
    assert result.skipped[0].kind == "excluded_related"
    assert result.skipped[0].flags == ("full_siblings",)
    assert result.skipped[0].reason == "Excluded: full siblings"


def test_visible_ranked_hides_zeroes_and_applies_top_without_mutating_result():
    rows = [
        _row("M01", "M", A="A/a"),
        _row("M02", "M", A="a/a"),
        _row("F01", "F", A="A/a"),
        _row("F02", "F", A="A/A"),
    ]
    result = _run(rows, build_run_settings(targets={"A": "a/a"}), ["A"])
    original_ranked = result.ranked

    shown, zero_count, hidden_by_top = visible_ranked(result, top=1)

    assert len(result.ranked) == 4
    assert result.ranked is original_ranked
    assert len(shown) == 1
    assert shown[0].probability == Fraction(1, 2)
    assert zero_count == 2
    assert hidden_by_top == 1


def test_result_currentness_tracks_relevant_settings_but_not_top():
    rows = [_row("M01", "M", A="A/a"), _row("F01", "F", A="A/a")]
    settings = build_run_settings(targets={"A": "A/a"}, top=1)
    result = _run(rows, settings, ["A"])
    current_signature = result.signature

    assert result_is_current(result, current_signature)
    assert result_is_current(result, result.signature)
    changed_target_signature = run_signature(
        ("upload", "tiny"),
        build_run_settings(targets={"A": "a/a"}, top=1),
    )
    display_only_top_signature = run_signature(
        ("upload", "tiny"),
        build_run_settings(targets={"A": "A/a"}, top=20),
    )
    assert not result_is_current(result, changed_target_signature)
    assert result_is_current(result, display_only_top_signature)


def test_failed_result_is_never_current():
    rows = [_row("M01", "M", A="A/a"), _row("F01", "F", A="A/a")]
    result = _run(rows, build_run_settings(targets={"A": "A++"}), ["A"])
    assert not result_is_current(result, result.signature)


def test_result_view_is_frozen_and_source_rows_are_unchanged():
    rows = [
        _row("M01", "M", A="A/a"),
        _row("F01", "F", A="A/a"),
    ]
    rows_before = deepcopy(rows)
    settings = build_run_settings(targets={"A": "A/a"})
    fieldnames = [*BASE_FIELDS, "A"]
    upload_args = build_loader_args(
        fieldnames,
        rows,
        MAPPING,
        first_data_row=2,
    )
    loader_rows_before = deepcopy(upload_args["rows"])
    result = run_ranking(
        upload_args,
        settings,
        gene_pool=["A"],
        fieldnames=fieldnames,
        column_mapping=MAPPING,
        base_signature=("upload", "tiny"),
    )

    assert rows == rows_before
    assert upload_args["rows"] == loader_rows_before
    with pytest.raises(FrozenInstanceError):
        result.error_message = "changed"
    with pytest.raises(AttributeError):
        result.ranked.append(result.ranked[0])


def test_skipped_count_categories_use_core_names():
    rows = [
        _row("M01", "M", A=""),
        _row("F01", "F", A="A/a"),
    ]
    result = _run(rows, build_run_settings(targets={"A": "A/a"}), ["A"])

    assert skipped_counts(result) == {
        "missing_genotype": 1,
        "missing_phase": 0,
        "excluded_related": 0,
        "all_lethal": 0,
    }


def test_skipped_counts_total_every_kind_in_multi_kind_run():
    rows = [
        _row("M01", "M", A=""),
        _row("M02", "M", A="A/a", mother_id="DAM", father_id="SIRE_M"),
        _row("F01", "F", A="A/a", mother_id="DAM", father_id="SIRE_F"),
    ]
    result = _run(
        rows,
        build_run_settings(targets={"A": "A/a"}, exclude_related=True),
        ["A"],
    )

    counts = skipped_counts(result)
    assert {pair.kind for pair in result.skipped} == {
        "missing_genotype",
        "excluded_related",
    }
    assert sum(counts.values()) == len(result.skipped)
    assert counts == {
        "missing_genotype": 1,
        "missing_phase": 0,
        "excluded_related": 1,
        "all_lethal": 0,
    }


def test_skipped_counts_include_unrecognized_observed_kind():
    result = RunResult(
        ranked=(),
        skipped=(
            SkippedPairView("M01", "F01", "future_kind", "example", ()),
        ),
        excluded_animals=(),
        notes=(),
        pedigree_warnings=(),
        signature=("sig",),
        error_message=None,
    )

    assert skipped_counts(result) == {
        "missing_genotype": 0,
        "missing_phase": 0,
        "excluded_related": 0,
        "all_lethal": 0,
        "future_kind": 1,
    }


def test_all_lethal_pair_is_skipped_while_other_pairs_rank():
    rows = [
        _row("M01", "M", A="A/A"),
        _row("M02", "M", A="A/a"),
        _row("F01", "F", A="A/A"),
        _row("F02", "F", A="A/a"),
    ]
    result = _run(
        rows,
        build_run_settings(targets={"A": "A/a"}, lethals={"A": "A/A"}),
        ["A"],
    )

    assert result.error_message is None
    assert len(result.ranked) == 3
    all_lethal = [
        pair for pair in result.skipped if pair.kind == "all_lethal"
    ]
    assert len(all_lethal) == 1
    assert (all_lethal[0].male_id, all_lethal[0].female_id) == ("M01", "F01")
    assert all_lethal[0].reason == (
        "All offspring are lethal for gene 'A' from male M01 and female F01."
    )


def test_missing_phase_pair_keeps_flags_while_other_pairs_rank():
    rows = [
        _row(
            "M01",
            "M",
            A="A/a",
            B="B/b",
            AB_phase="",
            mother_id="DAM",
            father_id="SIRE_M",
        ),
        _row(
            "M02",
            "M",
            A="A/a",
            B="B/b",
            AB_phase="A|B/a|b",
        ),
        _row(
            "F01",
            "F",
            A="A/a",
            B="B/b",
            AB_phase="A|B/a|b",
            mother_id="DAM",
            father_id="SIRE_F",
        ),
        _row(
            "F02",
            "F",
            A="A/a",
            B="B/b",
            AB_phase="A|B/a|b",
        ),
    ]
    settings = build_run_settings(
        targets={"A": "A/a"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
    )
    result = _run(rows, settings, ["A", "B"], ["AB_phase"])

    assert result.error_message is None
    assert any(pair.male_id == "M02" for pair in result.ranked)
    missing_phase = [
        pair for pair in result.skipped if pair.kind == "missing_phase"
    ]
    assert len(missing_phase) == 2
    m01_f01 = next(
        pair
        for pair in missing_phase
        if (pair.male_id, pair.female_id) == ("M01", "F01")
    )
    assert m01_f01.reason == (
        "Missing phase for animal 'M01' for linked gene pair ('A', 'B'). "
        "Phase must be written in the configured gene order for this pair."
    )
    assert m01_f01.flags == ("half_siblings_mother",)


def test_missing_phase_column_is_a_configuration_error_without_ranking():
    rows = [
        _row("M01", "M", A="A/a", B="B/b"),
        _row("F01", "F", A="A/a", B="B/b"),
    ]
    settings = build_run_settings(
        targets={"A": "A/a"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
    )

    result = _run(rows, settings, ["A", "B"])

    assert result.error_message is not None
    assert "Phase column 'AB_phase'" in result.error_message
    assert "not present in the uploaded headers" in result.error_message
    assert result.ranked == ()
    assert result.skipped == ()
