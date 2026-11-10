"""Locate per-study GWAS results and build the meta-analysis commands.

The pipeline writes one `*_y_allresults.tsv` per study arm, under a path that
encodes ancestry, modality (array or WGS) and analysis:

    progression/<ancestry>/[WGS/]<analysis>/<ancestry>_<analysis>_CHROMS_1_22/
        gwas_results/<model>/<ancestry>_<arm>_y_allresults.tsv

Non-NBA cohorts are all European array data and sit one level deeper:

    progression/non_NBA/<cohort>/<analysis>/EUR_<analysis>_CHROMS_1_22/...
"""

from __future__ import annotations

import glob
from pathlib import Path

import gwas_config as cfg

# Top-level folders under the progression root that are not ancestries.
NOT_ANCESTRIES = {'meta', 'non_NBA'}

# Modalities, as the path fragment that follows the ancestry folder.
NBA, WGS = '', 'WGS/'


def collect_result_paths(analysis: str, ancestries=None, modalities=(NBA, WGS),
                         include_non_nba: bool = True, prog_dir=None,
                         verbose: bool = True) -> list[str]:
    """Every per-study results file that should go into one meta-analysis.

    `ancestries=None` takes every ancestry folder present. Non-NBA cohorts are
    only included when the array modality and EUR are both selected.
    """
    prog_dir = Path(prog_dir or cfg.prep_path())
    model_dir = cfg.MODEL_DIRS[analysis]
    paths: list[str] = []

    found = set()
    for anc_dir in sorted(prog_dir.iterdir()):
        anc = anc_dir.name
        if not anc_dir.is_dir() or anc in NOT_ANCESTRIES:
            continue
        found.add(anc)
        if ancestries is not None and anc not in ancestries:
            continue
        for mod in modalities:
            paths += sorted(glob.glob(
                f'{prog_dir}/{anc}/{mod}{analysis}/{anc}_{analysis}_CHROMS_1_22'
                f'/gwas_results/{model_dir}/{anc}_*_y_allresults.tsv'))

    if ancestries is not None:
        missing = [a for a in ancestries if a not in found]
        if missing:
            print(f'WARNING: no folder for {missing}; available: {sorted(found)}')

    non_nba_dir = prog_dir / 'non_NBA'
    if (include_non_nba and non_nba_dir.is_dir() and NBA in modalities
            and (ancestries is None or 'EUR' in ancestries)):
        for cohort_dir in sorted(non_nba_dir.iterdir()):
            if not cohort_dir.is_dir():
                continue
            paths += sorted(glob.glob(
                f'{cohort_dir}/{analysis}/EUR_{analysis}_CHROMS_1_22'
                f'/gwas_results/{model_dir}/EUR_*_y_allresults.tsv'))

    paths = sorted(set(paths))

    if verbose:
        who = 'all ancestries' if ancestries is None else ', '.join(ancestries)
        mods = ', '.join(m.rstrip('/') or 'NBA' for m in modalities)
        print(f'{analysis} [{who}] [{mods}]: {len(paths)} files')
        for p in paths:
            print(f'  {study_label(p):28} {p}')

    return paths


def study_label(path) -> str:
    """Unique label for one results file.

    Basenames alone collide: WGS files repeat the array names, and every
    non-NBA CS_disc file is called EUR_case_control whatever the cohort.
    """
    p = Path(path)
    name = p.name.replace('_y_allresults.tsv', '')
    parts = p.parts
    if 'non_NBA' in parts:
        cohort = parts[parts.index('non_NBA') + 1]
        if cohort not in name:
            name = f'{name}_{cohort}'
    if 'WGS' in parts:
        name += ' (WGS)'
    return name


def labelled_paths(paths) -> dict[str, str]:
    """{study label: path}, ready to hand to the forest-plot helpers."""
    return {study_label(p): p for p in paths}


def plink_meta_command(paths, out_prefix) -> str:
    """PLINK 1.9 fixed- and random-effects meta-analysis over `paths`."""
    return (
        f'plink --meta-analysis {" ".join(map(str, paths))} + qt'
        ' --meta-analysis-snp-field ID'
        ' --meta-analysis-bp-field POS'
        ' --meta-analysis-a1-field A1'
        ' --meta-analysis-a2-field REF'
        ' --meta-analysis-chr-field "#CHROM"'
        ' --meta-analysis-se-field "SE"'
        f' --out {out_prefix}'
    )


def metal_command(paths, outdir, prefix, pipeline_dir=None, profile='gcb_gcp') -> str:
    """`run_metal.nf` invocation for the same set of per-study results."""
    pipeline_dir = Path(pipeline_dir or cfg.PIPELINE_DIR)
    return (
        f'nextflow run {pipeline_dir}/run_metal.nf -profile {profile}'
        f' --metal_input "{",".join(map(str, paths))}"'
        f' --metal_outdir "{outdir}"'
        f' --metal_prefix "{prefix}"'
    )
