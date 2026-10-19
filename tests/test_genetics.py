"""Tests for single-gene genetics, checked against textbook Mendelian ratios."""
from fractions import Fraction

import pytest

from colony_breeder.core import (
    AllOffspringLethalError,
    format_genotype,
    gametes,
    offspring_distribution,
    parse_genotype,
    probability_of_target,
)

P = parse_genotype  # short alias so expected values stay readable
F = Fraction


# ---------- parsing ----------

def test_genotype_order_does_not_matter():
    assert P("aA") == P("Aa")
    assert P("+/Cre") == P("Cre/+")


def test_both_styles_are_understood():
    assert P("Aa") == ("A", "a")
    assert P("Cre/+") == ("+", "Cre")
    assert P("ob/ob") == ("ob", "ob")
    assert P("Ay/a") == ("Ay", "a")


def test_surrounding_whitespace_is_ignored():
    assert P("  Aa  ") == P("Aa")
    assert P(" Cre / + ") == P("Cre/+")


def test_alleles_are_case_sensitive():
    assert P("AA") != P("aa")
    assert P("Aa") != P("aa")


def test_bad_genotypes_raise_a_helpful_error():
    for bad in ["", "A", "AaB", "A/", "/a", "A/a/b", "Cre+"]:
        with pytest.raises(ValueError, match="Cannot read genotype"):
            P(bad)


def test_format_genotype_uses_slash_style():
    assert format_genotype(P("aA")) == "A/a"
    assert format_genotype(P("+/Cre")) == "+/Cre"


# ---------- gametes ----------

def test_homozygote_makes_one_kind_of_gamete():
    assert gametes(P("AA")) == {"A": F(1)}


def test_heterozygote_makes_two_kinds_half_each():
    assert gametes(P("Aa")) == {"A": F(1, 2), "a": F(1, 2)}
    assert gametes(P("Cre/+")) == {"Cre": F(1, 2), "+": F(1, 2)}


# ---------- textbook crosses ----------

def test_Aa_x_Aa_is_1_2_1():
    dist = offspring_distribution("Aa", "Aa")
    assert dist == {P("AA"): F(1, 4), P("Aa"): F(1, 2), P("aa"): F(1, 4)}


def test_Aa_x_aa_is_half_half():
    dist = offspring_distribution("Aa", "aa")
    assert dist == {P("Aa"): F(1, 2), P("aa"): F(1, 2)}


def test_AA_x_aa_gives_all_Aa():
    assert offspring_distribution("AA", "aa") == {P("Aa"): F(1)}


def test_AA_x_AA_gives_all_AA():
    assert offspring_distribution("AA", "AA") == {P("AA"): F(1)}


def test_AA_x_Aa_is_half_half():
    dist = offspring_distribution("AA", "Aa")
    assert dist == {P("AA"): F(1, 2), P("Aa"): F(1, 2)}


def test_aa_x_aa_gives_all_aa():
    assert offspring_distribution("aa", "aa") == {P("aa"): F(1)}


def test_slash_style_crosses_work():
    dist = offspring_distribution("Cre/+", "+/+")
    assert dist == {P("Cre/+"): F(1, 2), P("+/+"): F(1, 2)}


def test_parent_order_does_not_change_the_result():
    assert offspring_distribution("aA", "Aa") == offspring_distribution("Aa", "aA")
    assert offspring_distribution("Aa", "aa") == offspring_distribution("aa", "Aa")


def test_probabilities_always_sum_to_exactly_one():
    genotypes = ["AA", "Aa", "aa"]
    for mum in genotypes:
        for dad in genotypes:
            assert sum(offspring_distribution(mum, dad).values()) == 1


def test_results_are_exact_fractions_not_floats():
    dist = offspring_distribution("Aa", "Aa", lethal=["AA"])
    for probability in dist.values():
        assert isinstance(probability, Fraction)
    assert dist[P("Aa")] == F(2, 3)  # exactly 2/3, not 0.6666...


# ---------- lethal genotypes ----------

def test_lethal_homozygote_is_removed_and_rest_renormalised():
    dist = offspring_distribution("Ay/a", "Ay/a", lethal=["Ay/Ay"])
    assert dist == {P("Ay/a"): F(2, 3), P("a/a"): F(1, 3)}


def test_lethal_that_cannot_occur_changes_nothing():
    dist = offspring_distribution("Ay/a", "a/a", lethal=["Ay/Ay"])
    assert dist == {P("Ay/a"): F(1, 2), P("a/a"): F(1, 2)}


def test_lethal_genotype_may_be_written_in_either_order():
    # a heterozygous lethal, spelled both ways round
    assert offspring_distribution("Aa", "Aa", lethal=["aA"]) == offspring_distribution(
        "Aa", "Aa", lethal=["Aa"]
    )


def test_a_single_lethal_string_is_accepted():
    dist = offspring_distribution("Ay/a", "Ay/a", lethal="Ay/Ay")
    assert dist == {P("Ay/a"): F(2, 3), P("a/a"): F(1, 3)}


def test_if_everything_is_lethal_we_raise_an_error():
    with pytest.raises(AllOffspringLethalError, match="lethal"):
        offspring_distribution("Ay/Ay", "Ay/Ay", lethal=["Ay/Ay"])


# ---------- probability of a target ----------

def test_target_aa_from_Aa_x_Aa_is_one_quarter():
    assert probability_of_target("Aa", "Aa", "aa") == F(1, 4)


def test_target_may_be_several_genotypes():
    # "dominant phenotype" = AA or Aa = 1/4 + 1/2
    assert probability_of_target("Aa", "Aa", ["AA", "Aa"]) == F(3, 4)


def test_target_that_cannot_be_produced_is_zero():
    assert probability_of_target("AA", "AA", "aa") == F(0)


def test_target_probability_respects_lethal_renormalisation():
    assert probability_of_target("Ay/a", "Ay/a", "a/a", lethal=["Ay/Ay"]) == F(1, 3)


def test_lethal_target_has_probability_zero():
    assert probability_of_target("Ay/a", "Ay/a", "Ay/Ay", lethal=["Ay/Ay"]) == F(0)
