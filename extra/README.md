# Supplemental analyses

Side analyses and reference steps run alongside the main GWAS. Nothing here is
part of the pipeline in the repo root - these are kept so the work can be
found again, re-run, or cited.

| Notebook | What it does | Relation to the pipeline |
|---|---|---|
| `locuszoom_prep.ipynb` | Exports GWAS, PLINK `.meta` and METAL `.TBL` results into LocusZoom format, aligning the METAL effect allele to ALT | produces the inputs for LocusZoom, `locuscompare_setaxes.R`, FUMA and the collider-bias notebook |
| `fuma_input_prep.ipynb` | Converts a LocusZoom export into FUMA's `rsID A1 A2 BETA SE N P` table | after `locuszoom_prep` |
| `prs_association.ipynb` | Tests a PD risk PRS against progression, with the GWAS's own covariates and PCs | after the GWAS has run |
| `collider_bias_correction.ipynb` | Index-event bias correction: SlopeHunter, Dudbridge 2019, Cai 2022 (**R kernel**) | after `locuszoom_prep` |

## Setup

These notebooks reuse the pipeline's configuration module, so the same `.env`
applies - see [../.env.template](../.env.template):

```
set -a; source ../.env; set +a
jupyter lab
```

The `gwas` package is imported from the parent directory by default; set
`GWAS_SCRIPTS_DIR` if these notebooks are run from somewhere that does not
have that relative layout.

`gcs_io.py` reads objects straight out of the bucket with the storage client,
for the notebooks that run without a gcsfuse mount. It needs
`GOOGLE_CLOUD_PROJECT` and `GWAS_WORK_BUCKET`.

Additional variables used only here:

| Variable | Used by |
|---|---|
| `GWAS_PRS_SSCORE` | `prs_association.ipynb` - PLINK2 `.sscore` file to test |
| `SH_INCIDENCE_FILE`, `SH_PROGNOSIS_FILE` | `collider_bias_correction.ipynb` - the two GWAS inputs |
| `SH_OUT_PREFIX`, `SH_REF_BFILE` | `collider_bias_correction.ipynb` - output prefix, 1000 Genomes PLINK prefix for LD clumping |

## Notes

**Dated filenames.** The curated clinical datasets and some result files carry
a build date in the name. Those are set as variables at the top of the
notebook that reads them (`DATASETS`, `IN_PATH`), with `YYYY-MM-DD`
placeholders - fill in the version actually being used.

**R packages** for the collider-bias notebook are installed from GitHub; the
install commands are in its first markdown cell.
