"""Core genetics logic.

Rule for this package: standard library only, and no knowledge of the CLI or
web app. Interfaces call into this package, never the other way around.
"""
from colony_breeder.core.colony import (
    Animal,
    load_colony,
    load_colony_rows,
    pedigree_warnings,
)
from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    format_genotype,
    gametes,
    offspring_distribution,
    parse_genotype,
    probability_of_target,
    remove_lethal,
)
from colony_breeder.core.ranking import (
    PairRankingResult,
    RankedPair,
    SkippedPair,
    generate_pairs,
    rank_pairs,
)
from colony_breeder.core.relatedness import pair_flags

__all__ = [
    "Animal",
    "AllOffspringLethalError",
    "PairRankingResult",
    "RankedPair",
    "SkippedPair",
    "format_genotype",
    "gametes",
    "generate_pairs",
    "load_colony",
    "load_colony_rows",
    "offspring_distribution",
    "parse_genotype",
    "pair_flags",
    "pedigree_warnings",
    "probability_of_target",
    "rank_pairs",
    "remove_lethal",
]
