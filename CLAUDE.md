# colony-breeder

Lab colony tool: given a spreadsheet of animals and genotypes, rank male x female
breeding pairs by the probability of producing a target genotype.

## How to work with me (IMPORTANT)
I'm a beginner-to-intermediate CS student building this as a portfolio project for
SWE internships. I want to understand everything.
- Build in small steps, ONE step at a time, then STOP and wait for my review.
- Explain your decisions. Keep code readable over clever.
- After each step tell me: what changed and why, how to run the tests, and one
  thing to try on my own.
- Flag biological assumptions (autosomal vs X-linked, no sex differences in
  recombination, etc.) so I can check them.
- Never write the whole project in one go.

## Core design rules
- Python 3.9+. Core genetics is standard-library only (no pandas in core).
- Core genetics is fully separate from interfaces. CLI and web app both call it.
- All probabilities use fractions.Fraction (exact; 1/3 stays 1/3).
- Genotypes: two alleles, "Aa" or slash style ("Cre/+", "Ay/a"). Order never matters.
- Lethal genotypes are removed and the rest renormalised to sum to 1.
- Multiple unlinked genes: probabilities multiply.
- Flag (don't hide) related animals (shared mother or father), with an option to
  exclude them. Also note pairs already in the same cage.
- Blank/unknown genotype: skip that pair, never guess.
- Every genetics rule gets a pytest test checked against textbook ratios.
- Do not claim no similar tool exists. README needs an honest Related Work section
  (MousebreedeR, Breeding Scheme Architect, GenomeMixer, a Swiss multigeneration
  planner, Transnetyx Colony, Labii). My angle: flexible column mapping for any
  spreadsheet + user-defined symbols + linked-gene probabilities in a pair-ranking tool.

## Build order
1. Project setup (DONE)
2. Core genetics, one autosomal gene (DONE)
3. Colony loading from CSV, rank all male x female pairs, multiple unlinked genes,
   sibling/cage flags, pedigree warnings, with tests (DONE)
4. CLI (--target GENE=GENOTYPE, --lethal, --top, --exclude-related)  <-- NEXT
5. X-linked genes (males have one copy)
6. Linked genes: recombination frequency per gene pair, each animal's phase
   (HB/hb vs Hb/hB), parental gametes (1-r)/2, recombinant gametes r/2,
   tests vs known values (10 map units = 10% recombinant)
7. Streamlit app: upload .xlsx/.csv without reformatting; dropdowns to map ID, sex,
   cage, mother, father, each gene; per-gene settings table for symbol meanings and
   target genotypes; lethal/sibling/linkage options; downloadable ranked table; note
   that uploaded data is not stored; Google Sheets via download or published link
8. README (what it does, input format, examples, math, limitations/roadmap,
   Related Work) + deploy on Streamlit Community Cloud
9. Later: multi-generation planning, breeder age limits, inbreeding coefficient

## Decisions already made in step 2 (src/colony_breeder/core/genetics.py)
- A genotype is a sorted tuple of two allele strings, e.g. ("A", "a"); sorting
  makes "aA" == "Aa". Output is displayed in slash style ("A/a").
- Lethal renormalisation is among surviving pups (Ay/a x Ay/a with Ay/Ay lethal
  gives Ay/a = 2/3, a/a = 1/3).
- probability_of_target accepts one genotype or a collection of genotypes.
- Alleles are case-sensitive strings.

## Decisions made in step 3c
- `rank_pairs` multiplies per-gene probabilities (Fractions). ASSUMES target genes are unlinked; linked genes are step 6.
- `RankedPair.per_gene` holds each gene's probability. A lethal list applies only to its own gene; a lethal entry for an untargeted gene is validated but ignored.
- A blank genotype in ANY targeted gene skips the pair; blanks in untargeted genes do not.
- Performance idea for later: cache probability by (gene, male genotype, female genotype).

## Decisions made in step 3d
- Flags: full_siblings, half_siblings_mother, half_siblings_father,
  parent_offspring, same_cage. Blank IDs/cages never match.
- half_siblings_* means "at least half" if a parent ID is blank.
- Exclusion covers all related flags except same_cage (CLI: --exclude-related).
- Checks one generation only; grandparent/cousin links are not detected
  (limitation; see step 9 inbreeding coefficient).
- Missing genotype skip takes precedence; skipped pairs keep their flags.

## Decisions made in step 3e
- pedigree_warnings returns (severity, message) entries; warnings identify
  detectable pedigree inconsistencies and info identifies parent IDs absent
  from the colony.
- Pedigree warnings are separate from pair ranking and do not change flags or
  probabilities.
- Parent loops are checked only between two animals; cycles involving three
  or more animals are not detected.
- A missing parent ID is informational, not a warning, because retired
  breeders are often absent from the current colony.

## Biological assumptions
- Genes are treated as autosomal until step 5 adds X-linked genetics.
- Target genes are treated as unlinked until step 6 adds linkage.
- Every animal is assumed to be a fertile, available breeder. For example,
  real `ob/ob` mice are infertile; the sample data ignores fertility.

## Environment
- Windows, PowerShell, Python 3.14, virtualenv in .venv
  (activate: .venv\Scripts\Activate.ps1). Install: pip install -e ".[dev]"
- Run tests: pytest
- Layout: src/colony_breeder/core/ (genetics), tests/, data/sample_mice.csv (fake data)
- Real lab data must never be committed (see .gitignore).
- Always run tests with: .venv\Scripts\python -m pytest   (the assistant's terminal does not auto-activate the venv)