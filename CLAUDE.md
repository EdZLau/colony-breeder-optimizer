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
4. CLI (--target GENE=GENOTYPE, --lethal, --top, --exclude-related) (DONE)
5. X-linked genes (males have one copy) (DONE)
6. Linked genes: recombination frequency per gene pair, each animal's phase
   (HB/hb vs Hb/hB), parental gametes (1-r)/2, recombinant gametes r/2,
   tests vs known values (10 map units = 10% recombinant)
   (DONE)
7. Streamlit app:
   - 7a: Read CSV/XLSX uploads into validated text rows (DONE).
      - 7a-i: `load_colony_rows` supports uploaded rows and source row offsets (DONE).
   - 7b-1: Add pure upload mapping/loading helpers and opt-in unknown-sex loading
     (DONE).
   - 7b-2: Build the upload UI for sheet/header selection, preview, column/sex
     mapping, loading, and warnings (DONE).
      - 7b-2a: Add pure presentation logic in `app/streamlit_logic.py` (DONE).
      - 7b-2b: Add the Streamlit page in `app/streamlit_app.py` and configure
        `.streamlit/config.toml` (DONE).
   - 7c: Add target, lethal, relatedness, and linkage settings and run core
     ranking. <-- NEXT
   - 7d: Show ranked results and warnings; provide a downloadable table and a
     clear note that uploaded data is not stored.
   - 7e: Deploy on Streamlit Community Cloud and document Google Sheets input
     via download or published link.
   - Latest full-suite result: 374 tests passed.
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
- `rank_pairs` multiplies per-gene probabilities (Fractions) for unpaired
  autosomal genes only; linked pairs are different (see step 6).
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
- Autosomal and X-linked inheritance are supported. Configured autosomal linked
  pairs use step 6 linkage; other autosomal genes and X-linked genes are
  treated as independent groups.
- X-linked genes are treated as non-pseudoautosomal; X-inactivation and
  sex-specific viability apart from listed lethal genotypes are not modeled.
- Every animal is assumed to be a fertile, available breeder. For example,
  real `ob/ob` mice are infertile; the sample data ignores fertility.

## Decisions made in step 4
- CLI is in `colony_breeder.cli`; `main(argv=None)` parses and prints, then
  calls the core APIs without implementing genetics.
- `--target` is required and repeatable; repeated values for one gene are
  alternative target genotypes. Target and lethal genes are deduplicated when
  selecting genotype columns to load.
- `--top` must be at least 1. Zero-probability pairs are hidden and counted;
  skipped pairs are categorized and show five examples unless
  `--show-skipped` is supplied.
- CLI errors are concise stderr messages with a nonzero exit code; invalid
  text encoding advises saving the file as UTF-8 CSV.
- Pedigree warnings and results go to stdout; errors go to stderr. Exit codes
  are 0 for success, 1 for data errors, and 2 for command-line usage errors.
- Zero-probability pairs are hidden by default. An untargeted autosomal lethal
  setting is ignored with an explanatory note; untargeted X-linked lethals
  affect joint survival (see step 5).
- The CLI uses default column names only; configurable column mapping is
  still deferred for the CLI. Step 7b-1 adds mapping helpers for the upload
  path; its selected gene columns are inspected but not loaded or parsed.

## Decisions made in step 5
- X-linked genes are marked by gene name with `x_linked_genes` in
  `load_colony` and `rank_pairs`, or repeatable CLI `--xlinked GENE`; there is
  no per-row marker. Males use `A/Y`, females use `A/a`; male `/y` is normalized
  to `/Y`. Internal allele tuples are sorted (so `Y` may sort before lowercase
  alleles); use `format_x_linked_genotype` to display the Y last.
- X-linked offspring use a joint (sex, genotype) distribution with a 1:1 sex
  ratio at conception. Lethals are removed and renormalized once over the
  joint distribution, including lethals on untargeted X-linked genes.
- For each sex `s`, `weight_s = 1/2 × product(surviving mass)` over all
  X-linked genes. `target_s` uses target mass for targeted genes and surviving
  mass for untargeted genes. Overall X-linked probability is
  `(target_F + target_M) / (weight_F + weight_M)`.
- X-linked `per_gene` probabilities are descriptive marginals among survivors
  and must not be multiplied. X-linked males transmit their single X without
  recombination. Double-heterozygous females can recombine between X-linked
  genes, but the current model ignores that and treats multiple X-linked genes
  as unlinked (`r = 1/2`). This is an approximation, especially for tightly
  linked genes; X-linked linkage is deferred.
- Animals with blank or unrecognized sex are excluded from pairing and reported
  through `excluded_animals` and `notes`. With the default
  `allow_unknown_sex=False`, loading raises `ValueError` for unusable sex when
  no X-linked gene is marked; when an X-linked gene is marked, the loader
  accepts the row and stores the uppercased raw sex text. Only the upload path
  opts into `allow_unknown_sex=True`, which stores blank or unrecognized sex as
  `""`. In every case `rank_pairs` excludes animals whose sex is not M or F.
- The CLI warns when an unmarked gene has genotypes in `X/Y` form. Untargeted
  X-linked lethals affect joint survival and are not ignored; untargeted
  autosomal lethals remain ignored with a note.
- `data/sample_mice.csv` includes an `Xgene` column with fake demonstration
  genotypes.

## Decisions made in step 6
- Linkage core (`core/linkage.py`) handles one pair of autosomal genes using
  exact `Fraction` probabilities. `r` is a recombination fraction, not a map
  distance. It is supplied by the user as a `Fraction`, string, or integer;
  floats and booleans are rejected, and `r` must be in `[0, 1/2]`. The same
  `r` is currently used for both sexes, although female mice typically
  recombine more than males. Parental gametes each have probability `(1-r)/2`;
  recombinant gametes each have probability `r/2`. At short distances, 10 map
  units is approximated by `r = 1/10`; map distance can exceed the observed
  recombination fraction because double crossovers can restore the parental
  allele arrangement.
- Phase is stored separately from genotypes on `Animal`, keyed by the
  configured ordered pair `(gene1, gene2)`. In phase text `A|B/a|b`, `|`
  separates the alleles at gene 1 and gene 2 on each chromosome, and `/`
  separates the two chromosomes. Alleles may have multiple characters, such
  as `Cre|Ay/+|a`. The phase column must follow the configured gene order.
  `parse_phased_genotype` parses and canonicalizes chromosome order;
  `validate_phase_matches_genotypes` checks that the phase agrees with both
  per-gene genotype cells.
- Phase is required only when an animal is heterozygous at both linked genes;
  otherwise it is inferred. The loader validates a nonblank phase when both
  genotype cells are present; if either genotype is blank, the phase is stored
  raw and unvalidated. At ranking time, supplied phases are validated again for
  every usable-sex animal with both genotypes, not only animals reached in a
  particular cross. A malformed or mismatched phase raises `ValueError` and
  aborts the whole ranking call. A missing required phase key raises a
  configuration error; a present but blank phase skips the pair as
  `missing_phase`. Missing genotype takes precedence.
- `load_colony` accepts `linked_pairs=[(gene1, gene2, r, phase_column)]`;
  `rank_pairs` accepts `linked_pairs=[(gene1, gene2, r)]`. Reject a gene
  repeated in or across pairs, an X-linked or unknown gene, or an invalid
  recombination frequency. An active linked pair has either gene targeted or
  a nonempty lethal list on either gene. Both genotypes at both genes are
  required for both animals, including an untargeted lethal partner. The
  configured phase column must exist even for an inactive linked pair. Loader
  and ranking configurations must use the same ordered gene pair.
- An active pair is calculated jointly, including lethal filtering and one
  renormalization over surviving offspring. Its joint target probability is
  multiplied once with independent autosomal factors and the X-linked group
  probability. Linked genes are removed from the independent-autosomal loop.
  A lethal on an untargeted partner still affects joint survival. If no
  offspring survive, the pair is skipped as `all_lethal`. An untargeted linked
  pair with lethals contributes a factor of 1 if any offspring survive. This
  differs from step 3c: an untargeted lethal on an unpaired autosomal gene is
  validated but ignored.
- Distinct linked pairs and unpaired genes are treated as independent groups;
  linked pairs are assumed to be on different chromosomes from each other and
  from unpaired genes.
- `RankedPair.linked_pair_probabilities` reports each configured active pair's
  joint factor. Targeted linked genes in `per_gene` are descriptive marginal
  probabilities among survivors, not factors to multiply, just like the
  descriptive X-linked marginals.
- Linked pairs are autosomal only. X-linked genes remain in their separate
  joint calculation and are treated as unlinked from one another (`r = 1/2`).
  X-linked males have one X and do not recombine when transmitting it: fathers
  transmit their X to daughters and Y to sons. Double-heterozygous females can
  recombine between X-linked genes, but this model ignores that; X-linked
  linkage is deferred and the unlinked assumption is approximate, worst for
  tightly linked genes. No pseudoautosomal region, X-inactivation, or
  sex-specific viability apart from listed lethals is modeled.
- Full-suite result after 6b-2: 264 tests passed.

## Decisions made in step 6c
- Repeatable `--linked GENE1,GENE2,R,PHASE_COLUMN` configures an autosomal pair.
  Gene names cannot contain commas. `R` is passed to the core as text; both
  `"1/10"` and `"0.1"` represent the exact fraction `1/10`. CLI help describes
  R as a recombination fraction from 0 to 0.5, not map distance.
- The CLI uses each option's ordered gene pair for both loading and ranking.
  The loader validates phase text against genotype columns and reports phase
  errors before printing results. If both genes use identical allele symbols,
  a swapped gene order in phase text cannot be detected from the genotypes.
- Linked-pair genes are included among CSV genotype columns. An untargeted
  lethal setting is not described as ignored for genes in linked pairs, since
  a lethal on either linked gene can affect joint survival.
- Output includes each active pair's exact joint factor from
  `RankedPair.linked_pair_probabilities`. Active lethal-only pairs are shown
  with a note that their factor affects survival only. Linked `per_gene`
  values are labeled descriptive marginals and must not be multiplied.
- Invalid CLI syntax and recombination text are usage errors (exit 2).
  Loading and ranking data/configuration errors are concise stderr messages
  with exit 1; warnings and results are printed to stdout.
- `data/sample_linked_FAKE.csv` contains fake demonstration records for trying
  linked CLI options; no existing sample data was changed.
- One recombination fraction `r` is used for both sexes, although female mice
  typically recombine more.
- Treating 10 map units as `r = 1/10` is an approximation: double crossovers
  can make map distance exceed the observed recombination fraction.
- Linked pairs and unpaired genes are treated as independent groups.
- Latest full-suite result after 6c: 280 tests passed.

## Decisions made in step 7b-1
- `UploadLoadResult` carries the loaded `animals`, warnings grouped by key,
  `unknown_sex_count`, `unknown_sex_ids`, selected-gene `distinct_values`, and
  an optional `error_message`.
- `guess_column_mapping(fieldnames)` guesses colony fields by normalizing
  headers to lowercase alphanumeric characters and matching known aliases.
  Optional cage and parent fields default to `"(none)"`; required animal ID
  and sex fields without a match return `None`.
- `guess_sex_value_mapping(rows, sex_column)` guesses only `m`/`male` as M and
  `f`/`female` as F. Every other distinct value, including blank, defaults to
  unknown (`""`).
- `summarize_distinct_values(rows, columns, limit=50)` lists distinct values in
  first-seen order, displays blanks as `"(blank)"`, and appends `"+N more"`
  when distinct values exceed the limit.
- `build_loader_args(...)` validates that selections still exist in the
  current headers and that mapped fields use distinct columns. `"(none)"` for
  cage or parent fields injects a blank column under a reserved name; a
  collision with that reserved name is an error. It returns copied rows, the
  mapping, the `first_data_row` offset, `gene_columns=[]`, and
  `allow_unknown_sex=True`.
- `load_from_upload(...)` validates selections and sex mappings, loads colony
  structure, collects upload and pedigree warnings, counts unknown-sex animals
  and returns their IDs, and summarizes selected gene-column values. It returns
  validation errors in `error_message`.
- Only `m`/`male` and `f`/`female` are guessed as M and F. Everything else,
  including blank values, defaults to unknown. Every distinct sex-column value
  must have an explicit mapping. A mapping for a value no longer present in
  the upload is stale and raises an error; a current value without a mapping
  also raises an error.
- `load_colony_rows` defaults to `allow_unknown_sex=False`. With no X-linked
  gene marked, the default loader still rejects unusable sex. Only the upload
  helper opts into `allow_unknown_sex=True`, which stores blank or
  unrecognized sex as `""`.
- `build_loader_args` returns `gene_columns=[]` on purpose, and
  `load_from_upload` consequently loads structure only. Genotypes are **not
  loaded in 7b-1**: selected gene columns are summarized for review but not
  parsed. Step 7c must reload the colony with target, lethal, X-linked, and
  linked genes, plus phase columns, before ranking. The commit messages record
  no rationale for the empty gene list; the evidence is the tests only,
  including the assertion that loader arguments contain no gene columns and
  the fake-CSV test that loads junk genotype text such as `"not-a-genotype"` in
  a selected gene column without error.
- Animal IDs are trimmed but duplicate checks are case-sensitive: `"m01"` and
  `"M01"` are different animals. This is a limitation if a lab expects
  case-insensitive uniqueness.
- `first_data_row` is the 1-based source row of the first record after the header (`header_row + 1`), and loader error row numbers count from it. Interior blank rows are kept as blank records in both the XLSX and CSV readers, and the loader skips them only after advancing its row counter, so later row numbers match the source. Trailing blank rows are trimmed. Tests: the XLSX and CSV interior-blank-row tests in `tests/test_upload.py` assert `Row 5 is missing an animal_id`. Not tested: for CSV, rows are counted as records, so a quoted cell containing a line break would make the row number differ from the physical line number. XLSX date/time warnings use worksheet row numbers.
- `ui_helpers.py` and `upload.py` do not import Streamlit; their upload and
  mapping logic can be tested independently of the UI.
- Upload warnings with structured column metadata are grouped by column name;
  warnings without a column use the `"Upload"` key. These keys can collide with
  app-generated groups named `"Pedigree"`, `"Cage"`, `"Relatedness"`, or
  `"Upload"` if an upload column has the same name.

## Decisions made in step 7b-2
- `default_column_selections(fieldnames, previous=None)` uses
  `guess_column_mapping`, keeping valid previous choices when possible.
  `"(none)"` is kept only for the optional `cage`, `mother_id`, and `father_id`
  fields; required fields must use real columns.
- `default_sex_value_mapping(rows, sex_column, previous=None)` preserves choices
  for observed raw sex values and guesses values that are newly observed.
  Mapping keys are raw strings, with `""` for blank cells. `sex_value_label`
  renders blank as `"(blank)"` for display only; that label is never a key.
- `sex_mapping_check(rows, sex_column, limit=20)` counts distinct values and
  marks the mapping table undisplayable when there are more than 20.
- `gene_column_options(fieldnames, column_mapping)` excludes columns mapped to
  colony fields. `drop_stale_genes(selected, options)` removes selections that
  are no longer available.
- `input_signature(...)` returns a hashable signature covering the file name,
  SHA-256 content hash (computed in the page), selected sheet, header row,
  column mapping, sex mapping, and gene selection. A loaded result is shown
  only while its saved signature matches the current inputs.
- `preview_rows(upload, limit=10)` returns a deep copy of the preview rows.
  `format_load_result(result)` returns a frozen display view containing counts,
  unknown-sex IDs, grouped warnings, errors, and selected-gene value summaries.
  `colony_sex_warnings(animals)` warns when no male or no female is loaded.
- The page reads uploads in memory and does not write temporary files, cache
  upload contents, or log row data. The visible note asks users with sensitive
  data to run the app locally.
- The page loads colony structure only. It does not load or validate genotypes
  or rank pairs, and it states that genotypes are NOT loaded and ranking is
  unavailable until step 7c.
- `.streamlit/config.toml` sets `server.maxUploadSize = 10` MB and
  `browser.gatherUsageStats = false`.
- Streamlit is optional and is listed in the `[app]` extra. The Streamlit page
  test uses `pytest.importorskip("streamlit")`, so the core test suite can run
  when the app dependency is not installed.
- Full-suite result after step 7b-2: 374 tests passed.

## Environment
- Windows, PowerShell, Python 3.14, virtualenv in .venv
  (activate: .venv\Scripts\Activate.ps1). Install: pip install -e ".[dev]"
- Run tests: pytest
- Layout: src/colony_breeder/core/ (genetics), tests/, data/sample_mice.csv (fake data)
- Real lab data must never be committed (see .gitignore).
- Always run tests with: .venv\Scripts\python -m pytest   (the assistant's terminal does not auto-activate the venv)