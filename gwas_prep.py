"""Turn a curated clinical dataset into the files the LongGWAS pipeline reads.

The three input-prep notebooks (NBA, WGS, non-NBA) differ only in how samples
are selected and which genotype files they point at; everything downstream of
that selection - covariate formatting, phenotype coding, the survival fixes and
the PLINK subsetting script - is shared and lives here.

Expected columns in the curated dataset (after `id` is renamed to GP2ID):
    GP2ID, cohort, sex, study_days, event, HY,
    age_at_baseline, age_at_onset, disease_duration,
    c_age, c_age_sq, c_disease_duration, c_dd_sq
plus study_type / GP2_phenotype, merged in from the release master key.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

# Clinical covariates per analysis, in the order they are written out.
# Cross-sectional models use the centred age / disease-duration terms; the
# survival models condition on time, so they use the raw baseline values.
CS_COVARS = ['SEX', 'PHENO', 'age_at_baseline', 'age_at_onset',
             'c_age', 'c_age_sq', 'disease_duration', 'c_disease_duration', 'c_dd_sq',
             'study_type', 'cohort']
SURV_COVARS = ['SEX', 'PHENO', 'age_at_baseline', 'disease_duration',
               'study_type', 'cohort']


def check_required_fields(subset: pd.DataFrame, analysis: str,
                          id_col: str = 'GP2ID') -> pd.DataFrame:
    """Drop samples missing a field the model needs, reporting what went.

    Cross-sectional models need an age at baseline and an age at onset; every
    model needs its outcome (HY for continuous CS, `event` otherwise).
    """
    subset = subset.copy()

    if analysis in ('CS', 'CS_disc', 'LT'):
        age = subset.groupby([id_col, 'age_at_baseline', 'age_at_onset'],
                             dropna=False, as_index=False).size()
        no_baseline = age[age.age_at_baseline.isnull()]
        no_onset = age[age.age_at_onset.isnull()]
        for name, bad in (('age_at_baseline', no_baseline), ('age_at_onset', no_onset)):
            if len(bad):
                print(f'dropping {bad[id_col].nunique()} samples with no {name}')
                subset = subset[~subset[id_col].isin(bad[id_col])]

    outcome_col = 'HY' if analysis in ('CS', 'LT') else 'event'
    if outcome_col in subset.columns:
        missing = subset[subset[outcome_col].isnull()]
        if len(missing):
            print(f'dropping {missing[id_col].nunique()} samples with no {outcome_col}')
            subset = subset[~subset[id_col].isin(missing[id_col])]
    else:
        print(f'WARNING: no {outcome_col} column in the input data')

    return subset


def format_covariates(subset: pd.DataFrame, analysis: str,
                      id_cols=('#FID', 'IID')) -> pd.DataFrame:
    """PLINK-style covariate table for one analysis.

    `subset` must already carry the PLINK ID columns alongside the clinical
    fields. SEX and PHENO are recoded to PLINK's conventions, and every
    categorical value is stripped of characters the pipeline splits on.
    """
    covar = subset.rename(columns={'sex': 'SEX', 'GP2_phenotype': 'PHENO'})
    covar = covar[list(id_cols) + (SURV_COVARS if analysis == 'SURV' else CS_COVARS)].copy()

    covar['PHENO'] = covar['PHENO'].replace({'PD': 1})
    covar['SEX'] = covar['SEX'].replace({'Male': 1, 'Female': 2})

    # Covariate names and levels are passed to the pipeline as whitespace- and
    # dash-separated strings, so neither may appear inside a value.
    covar = covar.replace({'Case(/Control)': 'case_control'})
    covar['cohort'] = covar['cohort'].str.replace('-', '_', regex=False)

    if analysis == 'LT':
        # The longitudinal model takes cohort as dummy columns rather than a
        # single categorical covariate.
        covar = pd.get_dummies(covar, columns=['cohort'], prefix='cohort')
        bool_cols = covar.select_dtypes(include=['bool']).columns
        covar[bool_cols] = covar[bool_cols].astype(int)

    return covar.drop_duplicates()


def format_phenotypes(subset: pd.DataFrame, analysis: str,
                      id_cols=('IID',)) -> pd.DataFrame:
    """Phenotype table for one analysis: ID columns, `y`, and `study_days`."""
    ids = list(id_cols)

    if analysis == 'CS':
        # Continuous HY at a single visit.
        outcome = subset[ids + ['HY']].rename(columns={'HY': 'y'})
        outcome['study_days'] = 0
    elif analysis == 'CS_disc':
        # Reaching HY >= 3 as a case/control outcome; PLINK wants 1/2, not 0/1.
        outcome = subset[ids + ['event']].rename(columns={'event': 'y'})
        outcome['study_days'] = 0
        outcome['y'] = outcome['y'].replace({0: 1, 1: 2})
    elif analysis == 'SURV':
        outcome = subset[ids + ['event', 'study_days']].rename(columns={'event': 'y'})
    elif analysis == 'LT':
        outcome = subset[ids + ['HY', 'study_days']].rename(columns={'HY': 'y'})
    else:
        raise ValueError(f'unknown analysis: {analysis}')

    outcome = outcome.drop_duplicates().dropna()

    repeats = outcome.groupby(ids, dropna=False, as_index=False).size()
    print(f'{len(outcome)} rows for {len(repeats)} samples '
          f'(mean {repeats["size"].mean():.2f} per sample)')
    print(outcome['y'].value_counts().to_string())
    return outcome


def drop_zero_time_events(covar: pd.DataFrame, pheno: pd.DataFrame,
                          id_col: str = 'IID'):
    """Survival only: an event recorded at t = 0 carries no follow-up, drop it."""
    keep = ~((pheno.study_days == 0) & (pheno.y == 1))
    print(f'dropping {(~keep).sum()} samples with an event at study_days == 0')
    pheno = pheno[keep]
    return covar[covar[id_col].isin(pheno[id_col])], pheno


def drop_small_cohorts(covar: pd.DataFrame, pheno: pd.DataFrame, min_n: int = 50,
                       id_col: str = 'IID'):
    """Drop cohorts too small to fit their own stratum in the survival model."""
    counts = covar.cohort.value_counts()
    small = counts[counts < min_n]
    if len(small):
        print(f'dropping {len(small)} cohorts with < {min_n} samples: '
              + ', '.join(f'{c} ({n})' for c, n in small.items()))
    covar = covar[~covar.cohort.isin(small.index)]
    return covar, pheno[pheno[id_col].isin(covar[id_col])]


def cohort_sample_ids(df: pd.DataFrame, cohort: str, id_map: pd.DataFrame | None = None,
                      double_id: bool = True, id_col: str = 'GP2ID') -> pd.DataFrame:
    """Map clinical IDs onto the sample IDs used in a non-NBA cohort's PLINK files.

    Most cohorts genotype under the clinical ID itself. Two need help:
      * OPDC clinical IDs look like "PDS/1234/5", the PLINK IDs like "12345";
      * PROBAND carries a separate map, written by the GenoTools notebook, with
        `patient_ID` and `IID` columns - pass it as `id_map`.

    `double_id` writes the "<IID>_<IID>" form that the normalised non-NBA psam
    files use; set it to False for cohorts whose psam keeps the bare IID.
    """
    # Drop any ID columns already present so a merge cannot produce IID_x/IID_y.
    out = df.drop(columns=[c for c in ('IID', '#IID') if c in df.columns])

    if id_map is not None:
        id_map = id_map.rename(columns={'patient_ID': id_col})
        out = out.merge(id_map[[id_col, 'IID']], on=id_col, how='inner')
    elif cohort == 'OPDC':
        out['IID'] = (out[id_col].str.split('PDS/').str[1]
                      .str.replace('/', '', regex=False))
    else:
        out['IID'] = out[id_col]

    out['#IID'] = out['IID'] + '_' + out['IID'] if double_id else out['IID']
    return out


def write_subset_script(script_path, pfile: str, keep_file: str, out_prefix: str,
                        dest: str | None = None, extra_args: str = '') -> Path:
    """Write a bash script that subsets one chromosome of a release pfile.

    The script takes the chromosome as its only argument, so `pfile` and
    `out_prefix` should contain the literal `${CHROM}` placeholder. When `dest`
    is given the output is copied there and the local copy removed, which keeps
    the VM's scratch disk from filling up over 22 chromosomes.
    """
    script_path = Path(script_path)
    script_path.parent.mkdir(parents=True, exist_ok=True)

    lines = ['#!/bin/bash', 'set -euo pipefail', '', 'CHROM=$1', '',
             f'plink2 --pfile {pfile} --keep {keep_file} {extra_args} '
             f'--make-pgen --out {out_prefix}'.replace('  ', ' ')]
    if dest:
        lines += [f'gsutil -m cp {out_prefix}.* {dest}', f'rm {out_prefix}.*']

    script_path.write_text('\n'.join(lines) + '\n')
    script_path.chmod(0o755)
    print(f'wrote {script_path}')
    return script_path


def keep_samples_across_analyses(base_dir, analyses, label: str,
                                 id_cols=('#FID', 'IID')) -> pd.DataFrame:
    """Union of the sample IDs used by several analyses of the same dataset.

    One genotype subset is shared by every analysis, so the keep-list has to
    cover all of them.
    """
    base_dir = Path(base_dir)
    frames = []
    for analysis in analyses:
        path = base_dir / analysis / f'hy_{analysis}_{label}_covariates.tsv'
        if not path.exists():
            print(f'missing file, skipping: {path}')
            continue
        frames.append(pd.read_csv(path, sep='\t', usecols=list(id_cols)))

    if not frames:
        raise FileNotFoundError(f'no covariate files found under {base_dir}')

    return (pd.concat(frames, ignore_index=True)
            .drop_duplicates(subset=list(id_cols))
            .reset_index(drop=True))
