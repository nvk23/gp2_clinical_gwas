"""Read files straight out of a GCS bucket, without a gcsfuse mount.

The main pipeline notebooks work off gcsfuse mounts; these supplemental
notebooks often run somewhere that has credentials but no mount, so they pull
objects through the storage client instead.

The project and bucket come from the environment (see the main pipeline's
.env.example):

    GOOGLE_CLOUD_PROJECT   project billed for the request
    GWAS_WORK_BUCKET       bucket holding the working data
"""

from __future__ import annotations

import os
from io import StringIO

import pandas as pd
from google.cloud import storage


def _project() -> str:
    project = os.environ.get('GWAS_GCS_PROJECT') or os.environ.get('GOOGLE_CLOUD_PROJECT', '')
    if not project:
        raise RuntimeError('set GOOGLE_CLOUD_PROJECT (or GWAS_GCS_PROJECT) '
                           'before reading from the bucket')
    return project


def get_bucket(bucket_name: str | None = None):
    """Handle on a requester-pays bucket, billed to this project."""
    bucket_name = bucket_name or os.environ.get('GWAS_WORK_BUCKET', '')
    if not bucket_name:
        raise RuntimeError('set GWAS_WORK_BUCKET or pass bucket_name')
    project = _project()
    client = storage.Client(project=project)
    return client.bucket(bucket_name, user_project=project)


def blob_as_csv(bucket, path: str, sep: str = r'\s+', header='infer') -> pd.DataFrame:
    """Download one object and parse it as a table."""
    blob = bucket.get_blob(path)
    if blob is None:
        raise FileNotFoundError(f'gs://{bucket.name}/{path}')
    text = blob.download_as_bytes().decode('utf-8')
    return pd.read_csv(StringIO(text), sep=sep, header=header)


def list_folders(bucket, prefix: str = '') -> list[str]:
    """Immediate subfolders of `prefix`, as full prefixes ending in '/'."""
    folders: set[str] = set()
    for page in bucket.list_blobs(prefix=prefix, delimiter='/').pages:
        folders.update(page.prefixes)
    return sorted(folders)
