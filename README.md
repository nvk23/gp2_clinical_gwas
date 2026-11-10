# GP2 clinical GWAS - analysis notebooks

Hoehn & Yahr progression GWAS: cross-sectional (reaching HY >= 3) and survival
(time to HY >= 3), run per ancestry and per cohort and combined by
meta-analysis.

## Order of operations

| Notebook | What it does | When to run it |
|---|---|---|
| `00_genotools_qc_imputation.ipynb` | Array QC, ancestry assignment, TOPMed imputation | **non-NBA cohorts only** - release cohorts arrive QCed and imputed |
| `01_prep_inputs_nba.ipynb` | Covariate / phenotype / genotype inputs for one ancestry, NBA array | once per ancestry x analysis |
| `02_prep_inputs_wgs.ipynb` | Same, for samples with WGS but no array data | once per ancestry x analysis |
| `03_prep_inputs_non_nba.ipynb` | Same, for one non-NBA cohort | once per cohort x analysis |
| `04_launch_gwas.ipynb` | Writes the params files and the Nextflow commands | after 01-03 have uploaded their files |
| `05_meta_analysis.ipynb` | Per-study QC, PLINK and METAL meta-analysis, Manhattan / QQ | after every dataset has run |
| `06_clean_results.ipynb` | Annotates results for GWASlab / LocusZoom, lead variants | after 05 |
| `07_forest_plots.ipynb` | Per-variant forest plots across studies | after 05 |
| `08_beta_beta_plots.ipynb` | Effect-size comparison between two analyses | after 06 |
| `locuscompare_setaxes.R` | LocusCompare figure for one locus, risk vs progression, on a shared axis | after 06, per locus |

Analyses are referred to throughout as `CS_disc` (cross-sectional, discrete
outcome, logistic regression) and `SURV` (survival, cox-proportional hazards). The pipeline also supports `CS` (continuous HY, linear regression)
and `LT` (longitudinal, linear mixed effects); the shared helpers handle all four, but final manuscript run
used only the first two.

## Setup

```
cp .env.template .env     # fill in buckets, prefixes and credentials
set -a; source .env; set +a
jupyter lab
```

No path is hard-coded: everything comes from the environment variables listed
in `.env.template`. `gwas_config.describe()` at the top of each notebook prints
what was resolved. Keep `.env` out of git.

## Shared modules

| Module | Contents |
|---|---|
| `gwas_config.py` | Environment-driven paths, bucket URLs, release layout, `shell_do` |
| `gwas_prep.py` | Covariate / phenotype formatting, survival fixes, PLINK subsetting |
| `gwas_launch.py` | Renders the params templates, builds the Nextflow commands |
| `gwas_meta.py` | Finds per-study results, builds the PLINK / METAL commands |
| `gwas_plots.py` | Manhattan, QQ, genomic inflation |
| `gwas_forest.py` | Forest plots, with or without a meta-analysis summary |
| `gwas_betabeta.py` | Beta-beta comparison and concordance statistics |

## Parameter templates

`templates/template_<modality>_<analysis>.yml` - one per modality
(`NBA`, `WGS`, `non_NBA`) and analysis (`CS_disc`, `SURV`), generalised from
the R12 run. Placeholders (`{{PREP_ROOT}}`, `{{ANCESTRY}}`, `{{COHORT}}`,
`{{GENO_ROOT}}`, `{{COVAR_NUMERIC}}`, `{{CHUNK_SIZE}}`, `{{ANALYSIS_NAME}}`)
are filled in by `gwas_launch.render_params()`; rendered files land in
`launch_ready/R<release>/`.

The model specification differs by analysis, not by dataset:

* **CS_disc** - centred age and disease-duration terms plus PC1-PC10, cohort
  as a categorical covariate, `study_arm_col: study_type`.
* **SURV** - raw baseline age and disease duration plus PC1-PC10, fitted per
  cohort (`study_arm_col: cohort`), `time_col: study_days`.

Override `covar_numeric` for a cohort too small to support ten PCs, and
`chunk_size` for a survival run that is being OOM-killed.

## locuscompare_setaxes.R

Run outside Jupyter, one locus at a time:

```
export LDLINK_TOKEN=...            # free: https://ldlink.nih.gov/?tab=apiaccess
export LC_RISK_FILE=...            # risk GWAS, LocusZoom format
export LC_PROGRESSION_FILE=...     # progression GWAS, LocusZoom format
export LC_OUT_DIR=...
Rscript locuscompare_setaxes.R
```

The locus (chromosome, position, alleles, lead SNP, window) and the plotting
options are set in the block at the top of the script. Inputs are the
bgzipped LocusZoom-format files written by
`extra/locuszoom_prep.ipynb`. All three panels share one
-log10(P) scale, which is the point of this version over stock
`locuscomparer`.
