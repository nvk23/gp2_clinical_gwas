"""Render LongGWAS parameter files and build the commands that launch them.

One template per (modality, analysis) lives in templates/; this module fills in
the placeholders and writes a ready-to-run params file under
launch_ready/R<release>/.
"""

from __future__ import annotations

import re
from pathlib import Path

import gwas_config as cfg

TEMPLATE_DIR = Path(__file__).resolve().parent / 'templates'
LAUNCH_DIR = Path(__file__).resolve().parent / 'launch_ready'

MODALITIES = ('NBA', 'WGS', 'non_NBA')

PCS = ' '.join(f'PC{i}' for i in range(1, 11))

# Default numeric covariates per analysis. Cohorts with few samples may need
# fewer PCs - pass covar_numeric= to override.
COVAR_NUMERIC = {
    'CS_disc': f'c_age c_age_sq c_disease_duration c_dd_sq {PCS}',
    'SURV': f'age_at_baseline disease_duration {PCS}',
}

# Variants per chunk. The Cox model holds the dosage matrix for a whole chunk
# in memory several times over, so large cohorts may need this lowered to stay
# inside the executor's memory limit (an OOM kill shows up as exit code 137).
DEFAULT_CHUNK_SIZE = 30000

_PLACEHOLDER = re.compile(r'\{\{(\w+)\}\}')


def render_params(analysis: str, modality: str, *, ancestry: str | None = None,
                  cohort: str | None = None, geno_root: str | None = None,
                  covar_numeric: str | None = None,
                  chunk_size: int = DEFAULT_CHUNK_SIZE,
                  analysis_name: str | None = None,
                  prep_root: str | None = None,
                  out_path=None, template_dir=None) -> Path:
    """Write one params YAML and return its path.

    `analysis` is CS_disc or SURV; `modality` is NBA, WGS or non_NBA. Non-NBA
    runs need `cohort` and `geno_root` (the gs:// folder with that cohort's
    normalised pgen files); the others need `ancestry`.
    """
    if modality not in MODALITIES:
        raise ValueError(f'modality must be one of {MODALITIES}')
    if modality == 'non_NBA':
        if not cohort or not geno_root:
            raise ValueError('non-NBA runs need cohort= and geno_root=')
        ancestry = ancestry or 'EUR'   # every non-NBA cohort here is European
    elif not ancestry:
        raise ValueError('pass ancestry= for NBA and WGS runs')

    template = Path(template_dir or TEMPLATE_DIR) / f'template_{modality}_{analysis}.yml'
    text = template.read_text()

    values = {
        'PREP_ROOT': prep_root or cfg.prep_url(),
        'GENO_ROOT': geno_root or '',
        'ANCESTRY': ancestry,
        'COHORT': cohort or '',
        'COVAR_NUMERIC': covar_numeric or COVAR_NUMERIC[analysis],
        'CHUNK_SIZE': str(chunk_size),
        # The pipeline names its output folder <analysis_name>_<analysis>_CHROMS_1_22;
        # downstream notebooks expect the ancestry here, including for non-NBA runs.
        'ANALYSIS_NAME': analysis_name or f'{ancestry}_CHROMS_CHR1_22',
    }

    rendered = _PLACEHOLDER.sub(lambda m: values[m.group(1)], text)
    left = _PLACEHOLDER.findall(rendered)
    if left:
        raise KeyError(f'unfilled placeholders in {template.name}: {sorted(set(left))}')

    if out_path is None:
        label = cohort if modality == 'non_NBA' else \
            (f'{ancestry}_WGS' if modality == 'WGS' else ancestry)
        out_path = LAUNCH_DIR / f'R{cfg.RELEASE}' / f'{label}_{analysis}.yml'
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered)
    print(f'wrote {out_path}')
    return out_path


def batch_env(store_root: str, project_name: str) -> str:
    """The `export` block the gcb2 Nextflow profile reads.

    Values come from the environment, so no credentials are stored in the
    notebook or in git. Print it and paste it into the shell that runs
    Nextflow.
    """
    cfg.require('TOWER_ORG', 'TOWER_WORKSPACE_ID', 'TOWER_ACCESS_TOKEN',
                'GOOGLE_CLOUD_PROJECT', 'GOOGLE_CLOUD_REGION',
                'GOOGLE_SERVICE_ACCOUNT_EMAIL', 'REFERENCE_DIR')
    return '\n'.join([
        f"export STORE_ROOT='{store_root}'",
        f"export PROJECT_NAME='{project_name}'",
        f"export TOWER_ORG='{cfg.TOWER_ORG}'",
        'export TOWER_ACCESS_TOKEN="$TOWER_ACCESS_TOKEN"',
        f"export REFERENCE_DIR='{cfg.REFERENCE_DIR}'",
        f"export GOOGLE_CLOUD_PROJECT='{cfg.GOOGLE_CLOUD_PROJECT}'",
        f"export GOOGLE_CLOUD_REGION='{cfg.GOOGLE_CLOUD_REGION}'",
        f"export GOOGLE_SERVICE_ACCOUNT_EMAIL='{cfg.GOOGLE_SERVICE_ACCOUNT_EMAIL}'",
        f"export TOWER_WORKSPACE_ID='{cfg.TOWER_WORKSPACE_ID}'",
    ])


def nextflow_command(params_file, pipeline_dir=None, profile='gcb2',
                     with_tower=True, resume=True) -> str:
    """`nextflow run` command for one params file."""
    parts = [f'nextflow run {Path(pipeline_dir or cfg.PIPELINE_DIR)}',
             f'-params-file {params_file}']
    if with_tower:
        parts.append('-with-tower')
    if resume:
        parts.append('-resume')
    parts.append(f'-profile {profile}')
    return ' \\\n    '.join(parts)
