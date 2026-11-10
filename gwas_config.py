"""Shared configuration for the GP2 clinical (HY progression) GWAS notebooks.

Every machine- and bucket-specific value is read from the environment so that
nothing user-specific ends up in the repository. Copy `.env.example`, fill it
in, and `source` it (or export the variables another way) before starting
Jupyter:

    set -a; source .env; set +a
    jupyter lab

`describe()` prints the resolved configuration, and `require()` fails early
with a readable message when a variable the notebook needs is still empty.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _env(name: str, default: str = '') -> str:
    return os.environ.get(name, default).strip()


# --- Local layout ----------------------------------------------------------
# Root of the analysis VM's scratch space; the gcsfuse mounts live underneath.
WORK_DIR = Path(_env('GWAS_WORK_DIR', str(Path.home() / 'gwas'))).expanduser()

TIER2_MOUNT = WORK_DIR / 'gp2tier2'   # read-only release data
WORK_MOUNT = WORK_DIR / 'working'     # read-write working bucket

# Fast local scratch for intermediate genotype files (deleted when the VM stops)
TMP_DIR = Path(_env('GWAS_TMP_DIR', '/dev/shm'))

# Reference files (hg38 fasta, liftOver chains, 1000G panel, TOPMed panel)
REF_DIR = Path(_env('GWAS_REF_DIR', str(WORK_DIR / 'ref_files')))

# Cloned tooling used by the non-NBA QC notebook
GENOTOOLS_DIR = Path(_env('GWAS_GENOTOOLS_DIR', str(WORK_DIR / 'GenoTools')))
QC_CODE_DIR = Path(_env('GWAS_QC_CODE_DIR', str(WORK_DIR / 'genetics-qc')))


# --- Buckets ---------------------------------------------------------------
# Bucket names only, without the gs:// scheme.
TIER2_BUCKET = _env('GWAS_TIER2_BUCKET')
WORK_BUCKET = _env('GWAS_WORK_BUCKET')

# Folder inside the working bucket that this project writes to.
PREFIX = _env('GWAS_PREFIX')

# Folder (inside WORK_MOUNT) holding the curated clinical datasets produced by
# the clinical data team. Relative to the mount, e.g. "<team>/Progression/datasets".
CLINICAL_DIR = _env('GWAS_CLINICAL_DIR')


# --- Release ---------------------------------------------------------------
RELEASE = int(_env('GWAS_RELEASE', '12'))

# The curated dataset filenames carry a build date, so they are configurable.
CLINICAL_FILES = {
    ('NBA', 'CS_disc'): _env('GWAS_CLIN_NBA_CS', 'clin_combined_cs_dataset_for_LR.csv'),
    ('NBA', 'SURV'): _env('GWAS_CLIN_NBA_SURV', 'clin_combined_surv_dataset_for_CPH.csv'),
    ('non_NBA', 'CS_disc'): _env('GWAS_CLIN_NON_NBA_CS', 'non_nba_clin_hy_cs_dataset_for_LR.csv'),
    ('non_NBA', 'SURV'): _env('GWAS_CLIN_NON_NBA_SURV', 'non_nba_clin_surv_dataset_for_CPH.csv'),
}

# Analyses the LongGWAS pipeline supports and the gwas_results subfolder each
# one writes to.
MODEL_DIRS = {'CS': 'glm', 'CS_disc': 'glm', 'LT': 'lmm_gallop', 'SURV': 'cph'}
ANALYSES = tuple(MODEL_DIRS)

# Autosomes, in the brace-expansion form the pipeline's `input` glob expects.
CHROM_GLOB = 'chr{' + ','.join(str(c) for c in range(1, 23)) + '}'


# --- Nextflow / Seqera (Tower) ---------------------------------------------
# Secrets stay in the environment; nothing here is ever written to a file.
TOWER_ORG = _env('GWAS_TOWER_ORG')
TOWER_WORKSPACE_ID = _env('GWAS_TOWER_WORKSPACE_ID')
TOWER_ACCESS_TOKEN = _env('TOWER_ACCESS_TOKEN')
GOOGLE_CLOUD_PROJECT = _env('GOOGLE_CLOUD_PROJECT')
GOOGLE_CLOUD_REGION = _env('GOOGLE_CLOUD_REGION')
GOOGLE_SERVICE_ACCOUNT_EMAIL = _env('GOOGLE_SERVICE_ACCOUNT_EMAIL')
PIPELINE_DIR = Path(_env('GWAS_PIPELINE_DIR', str(WORK_DIR / 'long-gwas-pipeline')))
REFERENCE_DIR = _env('GWAS_REFERENCE_DIR')  # gs:// URI used by the pipeline


def require(*names: str) -> None:
    """Raise if any of the named module-level settings is still empty."""
    missing = [n for n in names if not globals().get(n)]
    if missing:
        raise RuntimeError(
            'missing configuration: ' + ', '.join(missing)
            + ' - set the matching variables in your environment (see .env.example)'
        )


def describe() -> None:
    """Print the resolved configuration (no secrets)."""
    print(f'release        R{RELEASE}')
    print(f'work dir       {WORK_DIR}')
    print(f'tier2 mount    {TIER2_MOUNT}  <- gs://{TIER2_BUCKET or "?"}')
    print(f'working mount  {WORK_MOUNT}  <- gs://{WORK_BUCKET or "?"}')
    print(f'bucket prefix  {PREFIX or "?"}')
    print(f'clinical data  {CLINICAL_DIR or "?"}')
    print(f'scratch        {TMP_DIR}')
    print(f'tower token    {"set" if TOWER_ACCESS_TOKEN else "NOT SET"}')


# --- Paths -----------------------------------------------------------------
def mount_commands() -> list[str]:
    """gcsfuse commands that mount the release and working buckets."""
    require('TIER2_BUCKET', 'WORK_BUCKET')
    TIER2_MOUNT.mkdir(parents=True, exist_ok=True)
    WORK_MOUNT.mkdir(parents=True, exist_ok=True)
    return [
        f'gcsfuse --dir-mode 555 --file-mode 444 --implicit-dirs {TIER2_BUCKET} {TIER2_MOUNT}',
        f'gcsfuse --dir-mode 777 --file-mode 777 --implicit-dirs {WORK_BUCKET} {WORK_MOUNT}',
    ]


def gcs(*parts: str) -> str:
    """gs:// URL inside this project's folder of the working bucket."""
    require('WORK_BUCKET', 'PREFIX')
    return '/'.join([f'gs://{WORK_BUCKET}', PREFIX, *(str(p).strip('/') for p in parts)])


def as_gcs(path) -> str:
    """Translate a path under the local working mount into its gs:// URL."""
    require('WORK_BUCKET')
    rel = os.path.relpath(os.path.abspath(str(path)), WORK_MOUNT)
    if rel.startswith('..'):
        raise ValueError(f'{path} is not under {WORK_MOUNT}')
    return f'gs://{WORK_BUCKET}/{rel}'


def prep_parts(ancestry: str | None = None, cohort: str | None = None,
               wgs: bool = False) -> tuple[str, ...]:
    """Folder, relative to the progression root, holding one dataset's inputs.

    NBA ancestry -> ("EUR",); the same ancestry's WGS -> ("EUR", "WGS");
    a non-NBA cohort -> ("non_NBA", "PROBAND").
    """
    if cohort:
        return ('non_NBA', cohort)
    if not ancestry:
        raise ValueError('pass either ancestry= or cohort=')
    return (ancestry, 'WGS') if wgs else (ancestry,)


def prep_path(*parts: str) -> Path:
    """Local path under <working mount>/<prefix>/R<release>_prep/progression."""
    require('PREFIX')
    return WORK_MOUNT.joinpath(PREFIX, f'R{RELEASE}_prep', 'progression', *map(str, parts))


def prep_url(*parts: str) -> str:
    """gs:// equivalent of prep_path()."""
    return gcs(f'R{RELEASE}_prep', 'progression', *map(str, parts))


def clinical_path(analysis: str, source: str = 'NBA') -> Path:
    """Curated clinical dataset feeding one analysis ('CS_disc' or 'SURV')."""
    require('CLINICAL_DIR')
    key = (source, 'CS_disc' if analysis != 'SURV' else 'SURV')
    return WORK_MOUNT / CLINICAL_DIR / CLINICAL_FILES[key]


def master_key() -> Path:
    """Release master key (sample-level metadata and ancestry labels)."""
    return TIER2_MOUNT / f'release{RELEASE}' / 'clinical_data' / \
        f'master_key_release{RELEASE}_final_vwb.csv'


def imputed_pfile(ancestry: str, chrom: str | int = '${CHROM}') -> str:
    """PLINK2 pfile prefix for the imputed NBA genotypes of one ancestry."""
    return str(TIER2_MOUNT / f'release{RELEASE}' / 'imputed_genotypes' / ancestry /
               f'chr{chrom}_{ancestry}_release{RELEASE}_vwb')


def wgs_pfile(ancestry: str, chrom: str | int = '${CHROM}',
              caller: str = 'dragen_joint_calling') -> str:
    """PLINK2 pfile prefix for the joint-called WGS genotypes of one ancestry."""
    return str(TIER2_MOUNT / f'release{RELEASE}' / 'wgs' / caller / 'plink' / ancestry /
               f'chr{chrom}_{ancestry}_release{RELEASE}')


def covar_file(analysis: str, label: str) -> str:
    return f'hy_{analysis}_{label}_covariates.tsv'


def pheno_file(analysis: str, label: str) -> str:
    return f'hy_{analysis}_{label}_phenotypes.tsv'


# --- Shell -----------------------------------------------------------------
def shell_cmd(*parts) -> str:
    """Join command parts into a copy-pasteable, line-continued shell command."""
    return ' \\\n    '.join(str(p) for p in parts)


def shell_do(command: str, print_cmd: bool = True, log: bool = False,
             return_log: bool = False, err: bool = False):
    """Run `command` and optionally print or return its combined output."""
    if print_cmd:
        print(f'Executing: {" ".join(command.split())}', file=sys.stderr)

    res = subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')
