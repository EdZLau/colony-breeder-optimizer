# colony-breeder

Ranks male × female mouse breeding pairs by the exact probability that their pups
have a target genotype. Upload a colony spreadsheet, choose the genotype you want,
and get a ranked list of pairs.

**[Try the live app](https://colony-breeder-optimizer-m4y8lemtekyxfjzhzyp4su.streamlit.app/)**
(free hosting, so it may take a few seconds to wake up)

> Do not upload real lab data to the hosted app. Uploads are processed on the
> host's server. To use real data, run the app locally (instructions below).

## Try it in 2 minutes

1. Open the live app and upload `data/sample_mice.csv` from this repo (fake demo data).
2. Keep the guessed column and sex mappings. Under "Gene columns", pick `Cre`,
   then click **Load colony**.
3. Under "Breeding settings", set the target gene to `Cre` and the target
   genotype to `Cre/Cre`.
4. Click **Run ranking**.

You should see the top pairs at exactly 1/2, each flagged for how the animals are
related (parent-offspring, full siblings, half siblings).

## Status

Working: autosomal genes, X-linked genes, linked genes (with recombination
frequency), lethal genotypes, relatedness flags, a command-line interface, and a
Streamlit web app. The core genetics uses exact fractions, not floating point, and
is covered by 489 automated tests.

More sections (input format, the math, limitations, related work) are coming.
