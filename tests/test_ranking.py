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
    assert "F02" in skipped.reason
    assert "Cre" in skipped.reason


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
    assert all(animal_id in reason for animal_id in expected_animals)
    assert "Cre" in reason


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
    assert "M01" in reason
    assert "F02" in reason
    assert "Agouti" in reason
    assert "lethal" in reason


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
