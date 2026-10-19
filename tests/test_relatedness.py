from pathlib import Path

from colony_breeder.core.colony import Animal, load_colony
from colony_breeder.core.relatedness import pair_flags

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample_mice.csv"


def _animal(animal_id, sex, mother_id="", father_id="", cage=""):
    return Animal(
        animal_id=animal_id,
        sex=sex,
        cage=cage,
        mother_id=mother_id,
        father_id=father_id,
    )


def test_pair_flags_full_siblings_and_same_cage_are_sorted():
    male = _animal("M01", "M", mother_id="F00", father_id="M00", cage="C1")
    female = _animal("F01", "F", mother_id="F00", father_id="M00", cage="C1")

    assert pair_flags(male, female) == ["full_siblings", "same_cage"]


def test_pair_flags_distinguish_maternal_and_paternal_half_siblings():
    male = _animal("M01", "M", mother_id="F00", father_id="")
    female = _animal("F01", "F", mother_id="F00", father_id="M00")
    paternal_half_sibling = _animal("F02", "F", mother_id="F01", father_id="M00")

    assert pair_flags(male, female) == ["half_siblings_mother"]
    assert pair_flags(male, paternal_half_sibling) == []
    assert pair_flags(
        _animal("M02", "M", father_id="M00"),
        paternal_half_sibling,
    ) == ["half_siblings_father"]


def test_shared_mother_with_both_fathers_blank_is_maternal_half_siblings():
    male = _animal("M01", "M", mother_id="F00")
    female = _animal("F01", "F", mother_id="F00")

    assert pair_flags(male, female) == ["half_siblings_mother"]


def test_shared_mother_with_only_one_father_recorded_is_at_least_half_siblings():
    male = _animal("M01", "M", mother_id="F00", father_id="")
    female = _animal("F01", "F", mother_id="F00", father_id="M00")

    assert pair_flags(male, female) == ["half_siblings_mother"]


def test_blank_mothers_never_count_as_a_shared_mother():
    male = _animal("M01", "M", father_id="M00")
    female = _animal("F01", "F", father_id="M00")

    assert pair_flags(male, female) == ["half_siblings_father"]


def test_founders_with_blank_parent_ids_are_not_related():
    male = _animal("M01", "M")
    female = _animal("F01", "F")

    assert pair_flags(male, female) == []


def test_shared_parent_id_counts_even_when_parent_is_not_in_colony():
    male = _animal("M01", "M", father_id="unlisted-parent")
    female = _animal("F01", "F", father_id="unlisted-parent")

    assert pair_flags(male, female) == ["half_siblings_father"]


def test_pair_flags_parent_offspring_in_both_defined_directions():
    male_is_father = _animal("M01", "M")
    female_with_him_as_father = _animal("F01", "F", father_id="M01")
    female_is_mother = _animal("F02", "F")
    male_with_her_as_mother = _animal("M02", "M", mother_id="F02")

    assert pair_flags(male_is_father, female_with_him_as_father) == ["parent_offspring"]
    assert pair_flags(male_with_her_as_mother, female_is_mother) == ["parent_offspring"]


def test_whitespace_around_ids_is_ignored_for_parent_matching():
    male = _animal("M01 ", "M")
    female = _animal("F01", "F", father_id=" M01")

    assert pair_flags(male, female) == ["parent_offspring"]


def test_id_matching_is_case_sensitive():
    male = _animal("M01", "M", father_id="Parent")
    female = _animal("F01", "F", father_id="parent")

    assert pair_flags(male, female) == []


def test_cage_values_are_trimmed_before_matching():
    male = _animal("M01", "M", cage=" C1 ")
    female = _animal("F01", "F", cage="C1")

    assert pair_flags(male, female) == ["same_cage"]


def test_male_listed_as_another_animals_mother_does_not_crash():
    male = _animal("M01", "M", mother_id="M02")
    female = _animal("F01", "F")

    assert pair_flags(male, female) == []


def test_animal_listed_as_its_own_parent_does_not_crash_or_flag_unrelated_pair():
    male = _animal("M01", "M", mother_id="M01")
    female = _animal("F01", "F")

    assert pair_flags(male, female) == []


def test_same_id_recorded_as_mother_and_father_is_shared_in_both_roles():
    male = _animal("M01", "M", mother_id="P01", father_id="P01")
    female = _animal("F01", "F", mother_id="P01", father_id="P01")

    assert pair_flags(male, female) == ["full_siblings"]


def test_same_cage_flag_requires_a_nonblank_shared_cage():
    assert pair_flags(
        _animal("M01", "M", cage="C1"),
        _animal("F01", "F", cage="C1"),
    ) == ["same_cage"]
    assert pair_flags(_animal("M01", "M"), _animal("F01", "F")) == []


def test_spaces_only_cage_and_parent_ids_are_blank():
    male = _animal("M01", "M", mother_id="   ", father_id="  ", cage="   ")
    female = _animal("F01", "F", mother_id=" ", father_id=" ", cage=" ")

    assert pair_flags(male, female) == []


def test_sample_colony_relationship_and_cage_flags_match_expected_pairs():
    colony = load_colony(SAMPLE, gene_columns=["Cre"])
    males = [animal for animal in colony if animal.sex == "M"]
    females = [animal for animal in colony if animal.sex == "F"]
    actual = {
        flag: {
            (male.animal_id, female.animal_id)
            for male in males
            for female in females
            if flag in pair_flags(male, female)
        }
        for flag in (
            "full_siblings",
            "half_siblings_mother",
            "half_siblings_father",
            "parent_offspring",
            "same_cage",
        )
    }

    assert actual == {
        "full_siblings": {
            ("M03", "F03"),
            ("M03", "F04"),
            ("M04", "F05"),
        },
        "half_siblings_mother": {
            ("M05", "F03"),
            ("M05", "F04"),
        },
        "half_siblings_father": {
            ("M05", "F05"),
        },
        "parent_offspring": {
            ("M01", "F03"),
            ("M01", "F04"),
            ("M02", "F05"),
            ("M03", "F01"),
            ("M04", "F02"),
            ("M05", "F01"),
        },
        "same_cage": {
            ("M01", "F01"),
            ("M03", "F03"),
        },
    }
    m01_f01 = next(male for male in males if male.animal_id == "M01")
    f01 = next(female for female in females if female.animal_id == "F01")
    assert pair_flags(m01_f01, f01) == ["same_cage"]
    m03 = next(male for male in males if male.animal_id == "M03")
    f03 = next(female for female in females if female.animal_id == "F03")
    assert pair_flags(m03, f03) == ["full_siblings", "same_cage"]
