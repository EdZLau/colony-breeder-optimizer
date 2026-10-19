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


@pytest.mark.parametrize(
    ("male_id", "female_id", "target", "without_lethal", "with_lethal"),
    [
        ("M01", "F01", "a/Y", Fraction(1, 4), Fraction(0)),
        ("M02", "F01", "a/a", Fraction(1, 4), Fraction(1, 3)),
        ("M04", "F05", "a/a", Fraction(1, 2), Fraction(1)),
    ],
)
def test_sample_xgene_hand_calculated_probabilities(
    male_id, female_id, target, without_lethal, with_lethal
):
    colony = load_colony(
        SAMPLE,
        gene_columns=["Xgene"],
        x_linked_genes=["Xgene"],
    )

    plain = rank_pairs(
        colony,
        target={"Xgene": target},
        x_linked_genes=["Xgene"],
    )
    lethal = rank_pairs(
        colony,
        target={"Xgene": target},
        lethal={"Xgene": ["a/Y"]},
        x_linked_genes=["Xgene"],
    )

    def find_pair(result):
        return next(
            pair
            for pair in result.ranked
            if pair.male.animal_id == male_id
            and pair.female.animal_id == female_id
        )

    assert find_pair(plain).probability == without_lethal
    assert find_pair(lethal).probability == with_lethal


def test_x_linked_joint_target_and_lethal_probability_is_one_sixth():
    animals = [
        _animal("M01", "M", {"X1": "A/Y", "X2": "A/Y"}),
        _animal("F01", "F", {"X1": "A/a", "X2": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"X2": ["a/Y"]},
        lethal={"X1": ["a/Y"]},
        x_linked_genes=["X1", "X2"],
    )

    pair = result.ranked[0]
    assert pair.probability == Fraction(1, 6)
    assert pair.per_gene["X2"] == Fraction(1, 6)


def test_x_linked_target_probability_uses_joint_surviving_sex_branches():
    animals = [
        _animal("M01", "M", {"X1": "A/Y", "X2": "A/Y"}),
        _animal("F01", "F", {"X1": "A/a", "X2": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={
            "X1": ["A/A", "A/a"],
            "X2": ["A/a"],
        },
        lethal={"X1": ["A/Y", "a/Y"]},
        x_linked_genes=["X1", "X2"],
    )

    assert result.ranked[0].probability == Fraction(1, 2)
    assert result.ranked[0].per_gene == {
        "X1": Fraction(1),
        "X2": Fraction(1, 2),
    }


def test_untargeted_x_linked_lethal_shifts_another_genes_target_probability():
    animals = [
        _animal("M01", "M", {"X1": "A/Y", "X2": "A/Y"}),
        _animal("F01", "F", {"X1": "A/a", "X2": "A/a"}),
    ]

    without_lethal = rank_pairs(
        animals,
        target={"X2": "A/a"},
        x_linked_genes=["X1", "X2"],
    )
    with_untargeted_lethal = rank_pairs(
        animals,
        target={"X2": "A/a"},
        lethal={"X1": ["A/Y", "a/Y"]},
        x_linked_genes=["X1", "X2"],
    )

    assert without_lethal.ranked[0].probability == Fraction(1, 4)
    assert with_untargeted_lethal.ranked[0].probability == Fraction(1, 2)
    assert with_untargeted_lethal.ranked[0].per_gene["X2"] == Fraction(1, 2)


def test_mixed_sex_targets_across_x_linked_genes_have_zero_probability():
    animals = [
        _animal("M01", "M", {"X1": "a/Y", "X2": "a/Y"}),
        _animal("F01", "F", {"X1": "a/a", "X2": "a/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"X1": "a/a", "X2": "a/Y"},
        x_linked_genes=["X1", "X2"],
    )

    assert result.ranked[0].probability == Fraction(0)


def test_both_sex_targets_use_joint_survival_conditioning():
    animals = [
        _animal("M01", "M", {"X1": "A/Y", "X2": "A/Y"}),
        _animal("F01", "F", {"X1": "A/a", "X2": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={
            "X1": ["a/a", "a/Y"],
            "X2": ["a/a", "a/Y"],
        },
        lethal={"X1": ["A/A"], "X2": ["A/Y"]},
        x_linked_genes=["X1", "X2"],
    )

    # Joint survivors have mass 1/4 in each sex; the joint target mass is
    # 1/8 in sons, so its probability among survivors is 1/4.
    assert result.ranked[0].probability == Fraction(1, 4)


def test_x_linked_and_autosomal_targets_are_multiplied():
    animals = [
        _animal("M01", "M", {"X": "A/Y", "Cre": "Cre/+"}),
        _animal("F01", "F", {"X": "A/a", "Cre": "Cre/+"}),
    ]

    result = rank_pairs(
        animals,
        target={"X": "a/Y", "Cre": "Cre/Cre"},
        x_linked_genes=["X"],
    )

    assert result.ranked[0].per_gene == {
        "Cre": Fraction(1, 4),
        "X": Fraction(1, 4),
    }
    assert result.ranked[0].probability == Fraction(1, 16)


def test_per_gene_keys_follow_target_order_for_mixed_autosomal_and_x_linked_genes():
    animals = [
        _animal("M01", "M", {"X": "A/Y", "Cre": "Cre/+", "Ob": "ob/+"}),
        _animal("F01", "F", {"X": "A/a", "Cre": "Cre/+", "Ob": "ob/+"}),
    ]

    result = rank_pairs(
        animals,
        target={
            "X": "a/Y",
            "Ob": "ob/ob",
            "Cre": "Cre/Cre",
        },
        x_linked_genes=["X"],
    )

    assert list(result.ranked[0].per_gene) == ["X", "Ob", "Cre"]


@pytest.mark.parametrize(
    ("animal_id", "sex", "genotype"),
    [
        ("M01", "M", "A/a"),
        ("F01", "F", "A/Y"),
    ],
)
def test_rank_pairs_rejects_invalid_animal_x_linked_genotype(
    animal_id, sex, genotype
):
    animals = [
        _animal("M01", "M", {"X": "A/Y", "Cre": "Cre/+"}),
        _animal("F01", "F", {"X": "A/a", "Cre": "Cre/+"}),
    ]
    invalid_animal = next(animal for animal in animals if animal.animal_id == animal_id)
    invalid_animal.genotypes["X"] = genotype

    with pytest.raises(
        ValueError,
        match=rf"Invalid X-linked genotype for animal '{animal_id}' in gene 'X'",
    ):
        rank_pairs(
            animals,
            target={"Cre": "Cre/+"},
            x_linked_genes=["X"],
        )


def test_all_sons_lethal_but_daughters_survive_yields_zero_ranked_probability():
    animals = [
        _animal("M01", "M", {"X1": "A/Y", "X2": "A/Y"}),
        _animal("F01", "F", {"X1": "A/a", "X2": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"X2": "a/Y"},
        lethal={"X1": ["A/Y", "a/Y"]},
        x_linked_genes=["X1", "X2"],
    )

    assert len(result.ranked) == 1
    assert result.ranked[0].probability == Fraction(0)
    assert result.skipped == []


def test_all_x_linked_offspring_lethal_is_skipped_with_all_lethal_kind():
    animals = [
        _animal("M01", "M", {"X": "A/Y"}),
        _animal("F01", "F", {"X": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"X": "A/A"},
        lethal={"X": ["A/A", "A/a", "A/Y", "a/Y"]},
        x_linked_genes=["X"],
    )

    assert result.ranked == []
    assert len(result.skipped) == 1
    assert result.skipped[0].kind == "all_lethal"
    assert "X" in result.skipped[0].reason


def test_unusable_sex_animals_are_reported_but_never_paired():
    animals = [
        _animal("M01", "M", {"X": "A/Y", "Cre": "Cre/+"}),
        _animal("F01", "F", {"X": "A/a", "Cre": "+/+"}),
        _animal("U01", "?", {"X": "A/a", "Cre": "Cre/+"}),
    ]

    autosomal_result = rank_pairs(animals, target={"Cre": "Cre/+"})
    x_linked_result = rank_pairs(
        animals,
        target={"X": "A/a"},
        x_linked_genes=["X"],
    )

    for result in (autosomal_result, x_linked_result):
        assert len(result.ranked) == 1
        assert (result.ranked[0].male.animal_id, result.ranked[0].female.animal_id) == (
            "M01",
            "F01",
        )
        assert result.skipped == []
        assert len(result.excluded_animals) == 1
        assert result.excluded_animals[0].animal_id == "U01"
        assert result.excluded_animals[0].reason == (
            "sex is blank or unrecognized, so it can't be assigned a male or female role"
        )
        assert result.notes == [
            "Some animals were excluded because sex is blank or unrecognized, "
            "so they can't be assigned a male or female role."
        ]
    assert autosomal_result.notes == x_linked_result.notes


def test_missing_genotype_skip_precedes_relatedness_and_unknown_sex_is_reported():
    animals = [
        Animal("M01", "M", "", "F00", "M00", {"X": ""}),
        Animal("F01", "F", "", "F00", "M00", {"X": "A/a"}),
        _animal("U01", "", {"X": "A/a"}),
    ]

    result = rank_pairs(
        animals,
        target={"X": "A/a"},
        x_linked_genes=["X"],
        exclude_related=True,
    )

    assert result.ranked == []
    assert len(result.skipped) == 1
    assert result.skipped[0].kind == "missing_genotype"
    assert result.skipped[0].reason.startswith("Missing genotype")
    assert [item.animal_id for item in result.excluded_animals] == ["U01"]
