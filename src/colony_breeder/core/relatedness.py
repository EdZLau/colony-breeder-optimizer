"""Relationship and cage flags for a male-female pair."""

from typing import List

from colony_breeder.core.colony import Animal


def pair_flags(male: Animal, female: Animal) -> List[str]:
    """Return sorted relationship and same-cage flags for a breeding pair.

    A ``half_siblings_*`` flag means the pair shares that recorded parent:
    they are at least half siblings when the other parent's ID is blank.
    """
    male_mother = male.mother_id.strip()
    female_mother = female.mother_id.strip()
    male_father = male.father_id.strip()
    female_father = female.father_id.strip()
    male_id = male.animal_id.strip()
    female_id = female.animal_id.strip()
    male_cage = male.cage.strip()
    female_cage = female.cage.strip()

    # IDs are trimmed before comparison, but matching remains case-sensitive.
    shares_mother = bool(
        male_mother and female_mother and male_mother == female_mother
    )
    shares_father = bool(
        male_father and female_father and male_father == female_father
    )

    flags = []
    if male_cage and female_cage and male_cage == female_cage:
        flags.append("same_cage")
    if (
        (male_id and female_father and male_id == female_father)
        or (female_id and male_mother and female_id == male_mother)
    ):
        flags.append("parent_offspring")

    if shares_mother and shares_father:
        flags.append("full_siblings")
    else:
        if shares_mother:
            flags.append("half_siblings_mother")
        if shares_father:
            flags.append("half_siblings_father")

    return sorted(flags)
