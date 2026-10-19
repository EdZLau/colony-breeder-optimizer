"""Core genetics logic.

Rule for this package: standard library only, and no knowledge of the CLI or
web app. Interfaces call into this package, never the other way around.
"""
from colony_breeder.core.colony import Animal, load_colony
from colony_breeder.core.genetics import (
    format_genotype,
    gametes,
    offspring_distribution,
    parse_genotype,
    probability_of_target,
    remove_lethal,
)

__all__ = [
    "Animal",
    "format_genotype",
    "gametes",
    "load_colony",
    "offspring_distribution",
    "parse_genotype",
    "probability_of_target",
    "remove_lethal",
]
