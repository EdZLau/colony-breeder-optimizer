from fractions import Fraction
from pathlib import Path

import pytest

from colony_breeder.core.colony import Animal, load_colony
from colony_breeder.core.ranking import generate_pairs, rank_pairs

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample_mice.csv"


def _animal(animal_id, sex, genotypes):
    return Animal(
        animal_id=animal_id,
        sex=sex,
        cage="",
        mother_id="",
        father_id="",
        genotypes=genotypes,
    )


def test_generate_pairs_returns_every_male_by_female_combination():
    animals = [
        _animal("F02", "F", {"Cre": "+/+"}),
        _animal("M02", "M", {"Cre": "Cre/+"}),
        _animal("F01", "F", {"Cre": "+/+"}),
        _animal("M01", "M", {"Cre": "Cre/+"}),
    ]

    pairs = generate_pairs(animals)

    assert [(male.animal_id, female.animal_id) for male, female in pairs] == [
        ("M02", "F02"),
        ("M02", "F01"),
        ("M01", "F02"),
        ("M01", "F01"),
    ]


def test_rank_pairs_uses_exact_probabilities_and_sorts_ties_by_ids():
    animals = [
        _animal("F02", "F", {"Cre": "+/+"}),
        _animal("M02", "M", {"Cre": "Cre/+"}),
        _animal("F01", "F", {"Cre": "+/+"}),
        _animal("M01", "M", {"Cre": "Cre/+"}),
        _animal("M00", "M", {"Cre": "+/+"}),
    ]

    result = rank_pairs(animals, target={"Cre": "Cre/+"}, lethal={})

    assert [
        (pair.male.animal_id, pair.female.animal_id, pair.probability)
        for pair in result.ranked
    ] == [
        ("M01", "F01", Fraction(1, 2)),
        ("M01", "F02", Fraction(1, 2)),
        ("M02", "F01", Fraction(1, 2)),
        ("M02", "F02", Fraction(1, 2)),
        ("M00", "F01", Fraction(0)),
        ("M00", "F02", Fraction(0)),
    ]
    assert result.skipped == []


def test_rank_pairs_multiplies_probabilities_for_unlinked_target_genes():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Ob": "ob/+"}),
        _animal("F01", "F", {"Cre": "Cre/+", "Ob": "ob/+"}),
    ]

    result = rank_pairs(
        animals,
        target={"Cre": ["Cre/Cre"], "Ob": ["ob/ob"]},
    )

    assert len(result.ranked) == 1
    pair = result.ranked[0]
    assert pair.male.animal_id == "M01"
    assert pair.female.animal_id == "F01"
    assert pair.per_gene == {"Cre": Fraction(1, 4), "Ob": Fraction(1, 4)}
    assert pair.probability == Fraction(1, 16)
    assert pair.probability == pair.per_gene["Cre"] * pair.per_gene["Ob"]
    assert result.skipped == []


def test_rank_pairs_target_gene_order_does_not_change_probability():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Ob": "ob/+"}),
        _animal("F01", "F", {"Cre": "Cre/+", "Ob": "ob/+"}),
    ]

    first_order = rank_pairs(
        animals,
        target={"Cre": "Cre/Cre", "Ob": "ob/ob"},
    )
    reverse_order = rank_pairs(
        animals,
        target={"Ob": "ob/ob", "Cre": "Cre/Cre"},
    )

    assert first_order.ranked[0].probability == reverse_order.ranked[0].probability == Fraction(1, 16)
    assert first_order.ranked[0].per_gene == reverse_order.ranked[0].per_gene


def test_lethal_genotypes_affect_only_their_own_gene_probability():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Ob": "ob/+"}),
        _animal("F01", "F", {"Cre": "Cre/+", "Ob": "ob/+"}),
    ]

    result = rank_pairs(
        animals,
        target={"Cre": "Cre/Cre", "Ob": "ob/ob"},
        lethal={"Cre": ["Cre/+"]},
    )

    pair = result.ranked[0]
    assert pair.per_gene == {"Cre": Fraction(1, 2), "Ob": Fraction(1, 4)}
    assert pair.probability == Fraction(1, 8)
    assert pair.per_gene["Ob"] == Fraction(1, 4)


@pytest.mark.parametrize(
    ("target", "message"),
    [
        ({}, "at least one gene"),
        ({"Cre": []}, "genotype list.*Cre.*empty"),
    ],
)
def test_rank_pairs_rejects_empty_target_or_target_genotype_list(target, message):
    animals = [_animal("M01", "M", {"Cre": "Cre/+"})]

    with pytest.raises(ValueError, match=message):
        rank_pairs(animals, target=target)


def test_rank_pairs_skips_missing_target_genotype_with_animal_and_gene_reason():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Ob": "ob/+"}),
        _animal("F02", "F", {"Cre": "", "Ob": "ob/+"}),
    ]

    result = rank_pairs(animals, target={"Cre": "Cre/+"})

    assert result.ranked == []
    assert len(result.skipped) == 1
    skipped = result.skipped[0]
    assert (skipped.male.animal_id, skipped.female.animal_id) == ("M01", "F02")
    assert skipped.kind == "missing_genotype"
    assert "F02" in skipped.reason
    assert "Cre" in skipped.reason


def test_missing_genotype_skip_still_includes_relationship_and_cage_flags():
    animals = [
        Animal(
            animal_id="M01",
            sex="M",
            cage="C1",
            mother_id="F00",
            father_id="M00",
            genotypes={"Cre": ""},
        ),
        Animal(
            animal_id="F01",
            sex="F",
            cage="C1",
            mother_id="F00",
            father_id="M00",
            genotypes={"Cre": "Cre/+"},
        ),
    ]

    result = rank_pairs(animals, target={"Cre": "Cre/+"})

    assert result.ranked == []
    assert result.skipped[0].flags == ["full_siblings", "same_cage"]
    assert result.skipped[0].kind == "missing_genotype"


def test_missing_genotype_skip_takes_precedence_over_related_exclusion():
    animals = [
        Animal(
            animal_id="M01",
            sex="M",
            cage="C1",
            mother_id="F00",
            father_id="M00",
            genotypes={"Cre": ""},
        ),
        Animal(
            animal_id="F01",
            sex="F",
            cage="C1",
            mother_id="F00",
            father_id="M00",
            genotypes={"Cre": "Cre/+"},
        ),
    ]

    result = rank_pairs(
        animals,
        target={"Cre": "Cre/+"},
        exclude_related=True,
    )

    assert result.ranked == []
    assert len(result.skipped) == 1
    assert result.skipped[0].reason.startswith("Missing genotype")
    assert result.skipped[0].kind == "missing_genotype"
    assert result.skipped[0].flags == ["full_siblings", "same_cage"]


@pytest.mark.parametrize(
    ("male_genotype", "female_genotype", "expected_animals"),
    [
        ("", "Cre/+", ("M01",)),
        ("", "", ("M01", "F02")),
    ],
)
def test_rank_pairs_skips_blank_target_genotypes_and_names_animals_and_gene(
    male_genotype, female_genotype, expected_animals
):
    animals = [
        _animal("M01", "M", {"Cre": male_genotype}),
        _animal("F02", "F", {"Cre": female_genotype}),
    ]

    result = rank_pairs(animals, target={"Cre": "Cre/+"})

    assert result.ranked == []
    assert len(result.skipped) == 1
    reason = result.skipped[0].reason
    assert result.skipped[0].kind == "missing_genotype"
    assert all(animal_id in reason for animal_id in expected_animals)
    assert "Cre" in reason


def test_rank_pairs_names_every_blank_target_gene_and_animal_in_skip_reason():
    animals = [
        _animal("M01", "M", {"Cre": "", "Ob": "ob/+"}),
        _animal("F01", "F", {"Cre": "", "Ob": ""}),
    ]

    result = rank_pairs(
        animals,
        target={"Cre": "Cre/Cre", "Ob": "ob/ob"},
    )

    assert result.ranked == []
    assert len(result.skipped) == 1
    reason = result.skipped[0].reason
    assert "M01" in reason and "F01" in reason
    assert "Cre" in reason and "Ob" in reason
    assert "male M01 and female F01" in reason


@pytest.mark.parametrize(
    ("target", "lethal"),
    [
        ({"Cre": "Cre+"}, {}),
        ({"Cre": "Cre/+"}, {"Cre": ["Cre+"]}),
    ],
)
def test_rank_pairs_rejects_malformed_target_or_lethal_before_pair_ranking(target, lethal):
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+"}),
        _animal("F02", "F", {"Cre": "+/+"}),
    ]

    with pytest.raises(ValueError, match="Cannot read genotype"):
        rank_pairs(animals, target=target, lethal=lethal)


def test_blank_genotype_for_non_target_gene_does_not_skip_pair():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Ob": ""}),
        _animal("F02", "F", {"Cre": "+/+", "Ob": "ob/+"}),
    ]

    result = rank_pairs(animals, target={"Cre": "Cre/+"})

    assert len(result.ranked) == 1
    assert result.skipped == []
    assert result.ranked[0].probability == Fraction(1, 2)


def test_rank_pairs_applies_lethal_genotypes_for_the_target_gene():
    animals = [
        _animal("M01", "M", {"Agouti": "Ay/a"}),
        _animal("F02", "F", {"Agouti": "Ay/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"Agouti": "Ay/a"},
        lethal={"Agouti": ["Ay/Ay"]},
    )

    assert result.ranked[0].probability == Fraction(2, 3)


def test_lethal_setting_for_gene_outside_target_is_validated_but_does_not_affect_ranking():
    animals = [
        _animal("M01", "M", {"Cre": "Cre/+", "Agouti": "Ay/a"}),
        _animal("F01", "F", {"Cre": "+/+", "Agouti": "Ay/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"Cre": "Cre/+"},
        lethal={"Agouti": ["Ay/Ay"]},
    )

    assert len(result.ranked) == 1
    assert result.ranked[0].per_gene == {"Cre": Fraction(1, 2)}
    assert result.ranked[0].probability == Fraction(1, 2)
    assert result.skipped == []


def test_pair_with_all_offspring_lethal_is_skipped_with_animal_and_gene_reason():
    animals = [
        _animal("M01", "M", {"Agouti": "Ay/Ay"}),
        _animal("F02", "F", {"Agouti": "Ay/Ay"}),
    ]

    result = rank_pairs(
        animals,
        target={"Agouti": "Ay/a"},
        lethal={"Agouti": "Ay/Ay"},
    )

    assert result.ranked == []
    assert len(result.skipped) == 1
    reason = result.skipped[0].reason
    assert result.skipped[0].kind == "all_lethal"
    assert "M01" in reason
    assert "F02" in reason
    assert "Agouti" in reason
    assert "lethal" in reason


def test_exclude_related_skips_exactly_twelve_sample_pairs_with_flags():
    colony = load_colony(SAMPLE, gene_columns=["Cre"])

    result = rank_pairs(
        colony,
        target={"Cre": ["Cre/Cre"]},
        exclude_related=True,
    )

    excluded_pairs = {
        (pair.male.animal_id, pair.female.animal_id)
        for pair in result.skipped
        if pair.reason.startswith("Excluded:")
    }
    assert len(excluded_pairs) == 12
    assert all(
        pair.kind == "excluded_related"
        for pair in result.skipped
        if pair.reason.startswith("Excluded:")
    )
    assert all(
        pair.kind == "excluded_related"
        for pair in result.skipped
        if pair.reason.startswith("Excluded:")
    )
    assert excluded_pairs == {
        ("M03", "F03"),
        ("M03", "F04"),
        ("M04", "F05"),
        ("M05", "F03"),
        ("M05", "F04"),
        ("M05", "F05"),
        ("M01", "F03"),
        ("M01", "F04"),
        ("M02", "F05"),
        ("M03", "F01"),
        ("M04", "F02"),
        ("M05", "F01"),
    }
    assert all(pair.flags for pair in result.skipped if pair.reason.startswith("Excluded:"))
    excluded = [pair for pair in result.skipped if pair.reason.startswith("Excluded:")]
    assert sum("full_siblings" in pair.flags for pair in excluded) == 3
    assert sum("half_siblings_mother" in pair.flags for pair in excluded) == 2
    assert sum("half_siblings_father" in pair.flags for pair in excluded) == 1
    assert sum("parent_offspring" in pair.flags for pair in excluded) == 6
    assert (
        sum("full_siblings" in pair.flags for pair in excluded)
        + sum("half_siblings_mother" in pair.flags for pair in excluded)
        + sum("half_siblings_father" in pair.flags for pair in excluded)
        + sum("parent_offspring" in pair.flags for pair in excluded)
        == 12
    )
    assert {pair.reason for pair in result.skipped if pair.reason.startswith("Excluded:")} == {
        "Excluded: full siblings",
        "Excluded: maternal half siblings",
        "Excluded: paternal half siblings",
        "Excluded: parent-offspring",
    }

    same_cage_only_pair = next(
        pair
        for pair in result.ranked
        if pair.male.animal_id == "M01" and pair.female.animal_id == "F01"
    )
    assert same_cage_only_pair.flags == ["same_cage"]


@pytest.mark.parametrize(
    ("target", "lethal"),
    [
        ({"Missing": "Aa"}, {}),
        ({"Cre": "Cre/+"}, {"Missing": "aa"}),
    ],
)
def test_rank_pairs_rejects_gene_names_not_in_colony(target, lethal):
    animals = [_animal("M01", "M", {"Cre": "Cre/+"})]

    with pytest.raises(ValueError, match="not present in the colony.*Missing"):
        rank_pairs(animals, target=target, lethal=lethal)


def test_sample_colony_cre_target_has_expected_skips_and_top_pairs():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    result = rank_pairs(colony, target={"Cre": ["Cre/Cre"]})

    assert len(result.ranked) == 30
    assert len(result.skipped) == 5
    assert {
        animal_id
        for pair in result.skipped
        for animal_id in (pair.male.animal_id, pair.female.animal_id)
        if animal_id.startswith("F")
    } == {"F06"}
    assert [
        (pair.male.animal_id, pair.female.animal_id, pair.probability)
        for pair in result.ranked[:3]
    ] == [
        ("M01", "F04", Fraction(1, 2)),
        ("M03", "F04", Fraction(1, 2)),
        ("M05", "F04", Fraction(1, 2)),
    ]


def test_sample_colony_ob_target_skips_no_pairs_despite_blank_cre_genotype():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    result = rank_pairs(colony, target={"Ob": ["ob/ob"]})

    assert len(result.ranked) == 35
    assert result.skipped == []


def test_sample_colony_cre_and_ob_targets_have_expected_top_pairs():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    result = rank_pairs(
        colony,
        target={"Cre": ["Cre/Cre"], "Ob": ["ob/ob"]},
    )

    assert len(result.ranked) == 30
    assert len(result.skipped) == 5
    assert all(
        pair.male.animal_id == "F06" or pair.female.animal_id == "F06"
        for pair in result.skipped
    )
    assert [
        (pair.male.animal_id, pair.female.animal_id, pair.probability)
        for pair in result.ranked[:4]
    ] == [
        ("M03", "F04", Fraction(1, 4)),
        ("M01", "F04", Fraction(1, 8)),
        ("M03", "F01", Fraction(1, 8)),
        ("M01", "F01", Fraction(1, 16)),
    ]


def test_sample_colony_agouti_lethal_setting_adjusts_top_pair_probabilities():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    result = rank_pairs(
        colony,
        target={
            "Cre": ["Cre/Cre"],
            "Ob": ["ob/ob"],
            "Agouti": ["a/a"],
        },
        lethal={"Agouti": ["Ay/Ay"]},
    )

    assert [
        (pair.male.animal_id, pair.female.animal_id, pair.probability)
        for pair in result.ranked[:4]
    ] == [
        ("M03", "F04", Fraction(1, 12)),
        ("M01", "F04", Fraction(1, 24)),
        ("M03", "F01", Fraction(1, 24)),
        ("M01", "F01", Fraction(1, 48)),
    ]
    assert result.ranked[0].per_gene == {
        "Cre": Fraction(1, 2),
        "Ob": Fraction(1, 2),
        "Agouti": Fraction(1, 3),
    }


def test_sample_colony_ob_and_agouti_targets_skip_no_pairs():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    result = rank_pairs(
        colony,
        target={"Ob": ["ob/ob"], "Agouti": ["a/a"]},
    )

    assert len(result.ranked) == 35
    assert result.skipped == []
