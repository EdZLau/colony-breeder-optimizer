"""Hand-checked tests for linked autosomal gene calculations."""
from fractions import Fraction

import pytest

from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    parse_genotype,
    probability_of_target,
)
from colony_breeder.core.linkage import (
    linked_offspring_distribution,
    parse_phased_genotype,
    phase_for_unambiguous_genotypes,
    phase_is_required,
    probability_of_linked_target,
    validate_phase_matches_genotypes,
)

F = Fraction
P = parse_genotype


def test_parse_phased_genotype_canonicalizes_chromosome_order():
    assert parse_phased_genotype("A|B/a|b") == parse_phased_genotype("a|b/A|B")
    assert parse_phased_genotype("A|B/a|b") == (("A", "B"), ("a", "b"))


def test_parse_phased_genotype_supports_multicharacter_alleles():
    assert parse_phased_genotype("Cre|+/+|+") == (("+", "+"), ("Cre", "+"))


@pytest.mark.parametrize(
    "malformed",
    [
        "",
        "A|B",
        "A|B/a|b/c|d",
        "A/a|b",
        "A|B/a|",
        "A||B/a|b",
    ],
)
def test_parse_phased_genotype_rejects_malformed_text(malformed):
    with pytest.raises(ValueError, match="phased genotype"):
        parse_phased_genotype(malformed)


def test_validate_phase_checks_against_per_gene_genotypes():
    validate_phase_matches_genotypes("A|B/a|b", "A/a", "B/b")

    with pytest.raises(ValueError, match="does not agree"):
        validate_phase_matches_genotypes("A|b/a|B", "A/A", "B/b")


@pytest.mark.parametrize(
    ("gene1", "gene2", "required"),
    [
        ("A/a", "B/b", True),
        ("A/A", "B/b", False),
        ("A/a", "B/B", False),
        ("A/A", "B/B", False),
    ],
)
def test_phase_is_required_only_when_both_genes_are_heterozygous(
    gene1, gene2, required
):
    assert phase_is_required(gene1, gene2) is required


@pytest.mark.parametrize(
    ("frequency", "expected"),
    [
        (F(0), F(1, 2)),
        (F(1, 2), F(1, 4)),
        ("1/10", F(9, 20)),
        ("0.1", F(9, 20)),
    ],
)
def test_coupling_phase_test_cross_has_hand_calculated_targets(frequency, expected):
    assert probability_of_linked_target(
        "A|B/a|b",
        "a|b/a|b",
        frequency,
        target_gene1="A/a",
        target_gene2="B/b",
    ) == expected


def test_repulsion_phase_changes_the_result_at_one_tenth_recombination():
    assert probability_of_linked_target(
        "A|b/a|B",
        "a|b/a|b",
        "1/10",
        target_gene1="A/a",
        target_gene2="B/b",
    ) == F(1, 20)


def test_half_recombination_matches_product_of_single_gene_probabilities():
    linked = probability_of_linked_target(
        "A|B/a|b",
        "a|b/a|b",
        F(1, 2),
        target_gene1="A/a",
        target_gene2="B/b",
    )
    independent = probability_of_target("Aa", "aa", "Aa") * probability_of_target(
        "Bb", "bb", "Bb"
    )

    assert linked == independent == F(1, 4)


def test_r_zero_double_heterozygote_cross_has_three_joint_outcomes():
    assert linked_offspring_distribution(
        "A|B/a|b",
        "A|B/a|b",
        F(0),
    ) == {
        (P("A/A"), P("B/B")): F(1, 4),
        (P("A/a"), P("B/b")): F(1, 2),
        (P("a/a"), P("b/b")): F(1, 4),
    }


def test_lethal_at_one_linked_gene_renormalizes_joint_distribution():
    assert linked_offspring_distribution(
        "A|B/a|b",
        "A|B/a|b",
        F(0),
        lethal_gene1="a/a",
    ) == {
        (P("A/A"), P("B/B")): F(1, 3),
        (P("A/a"), P("B/b")): F(2, 3),
    }


def test_untargeted_partner_lethal_changes_target_probability():
    without_lethal = probability_of_linked_target(
        "A|B/a|b",
        "a|b/a|b",
        "1/10",
        target_gene1="A/a",
    )
    with_untargeted_partner_lethal = probability_of_linked_target(
        "A|B/a|b",
        "a|b/a|b",
        "1/10",
        target_gene1="A/a",
        lethal_gene2="B/b",
    )

    assert without_lethal == F(1, 2)
    assert with_untargeted_partner_lethal == F(1, 10)


def test_all_joint_offspring_lethal_raises_existing_error():
    with pytest.raises(AllOffspringLethalError, match="lethal"):
        linked_offspring_distribution(
            "A|B/a|b",
            "A|B/a|b",
            F(0),
            lethal_gene1=["A/A", "A/a", "a/a"],
        )


@pytest.mark.parametrize(
    "frequency",
    [
        F(-1, 10),
        F(3, 5),
        "not a fraction",
        "1/0",
        0.1,
        True,
        False,
    ],
)
def test_invalid_recombination_frequency_is_rejected(frequency):
    with pytest.raises(ValueError, match="recombination frequency"):
        linked_offspring_distribution("A|B/a|b", "a|b/a|b", frequency)


@pytest.mark.parametrize("frequency", [F(0), "1/10", F(1, 2)])
def test_homozygous_parent_cross_has_same_distribution_for_any_recombination(
    frequency,
):
    phase = phase_for_unambiguous_genotypes("A/A", "B/b")
    assert phase == parse_phased_genotype("A|B/A|b")
    assert linked_offspring_distribution(
        phase,
        "a|b/a|b",
        frequency,
    ) == {
        (P("A/a"), P("B/b")): F(1, 2),
        (P("A/a"), P("b/b")): F(1, 2),
    }


@pytest.mark.parametrize(
    ("gene1", "gene2", "expected"),
    [
        ("A/A", "B/B", parse_phased_genotype("A|B/A|B")),
        ("A/A", "B/b", parse_phased_genotype("A|B/A|b")),
        ("A/a", "B/B", parse_phased_genotype("A|B/a|B")),
    ],
)
def test_phase_for_unambiguous_genotypes_builds_phase(gene1, gene2, expected):
    assert phase_for_unambiguous_genotypes(gene1, gene2) == expected


def test_phase_for_unambiguous_genotypes_rejects_two_heterozygous_genes():
    with pytest.raises(ValueError, match="phase must be supplied"):
        phase_for_unambiguous_genotypes("A/a", "B/b")


@pytest.mark.parametrize("frequency", [F(0), F(1, 10), F(1, 2)])
def test_linked_distribution_sums_to_one_before_and_after_lethal_removal(
    frequency,
):
    distribution = linked_offspring_distribution(
        "A|B/a|b",
        "A|B/a|b",
        frequency,
    )
    surviving_distribution = linked_offspring_distribution(
        "A|B/a|b",
        "A|B/a|b",
        frequency,
        lethal_gene1="a/a",
    )

    assert sum(distribution.values(), F(0)) == F(1)
    assert sum(surviving_distribution.values(), F(0)) == F(1)


@pytest.mark.parametrize("frequency", [0, 1])
def test_integer_recombination_frequency_is_accepted_or_range_checked(frequency):
    if frequency == 0:
        assert probability_of_linked_target(
            "A|B/a|b",
            "a|b/a|b",
            frequency,
            target_gene1="A/a",
            target_gene2="B/b",
        ) == F(1, 2)
    else:
        with pytest.raises(ValueError, match="between 0 and 1/2"):
            linked_offspring_distribution("A|B/a|b", "a|b/a|b", frequency)
