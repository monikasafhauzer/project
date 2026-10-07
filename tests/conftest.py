import pytest


@pytest.fixture(autouse=True)
def isolated_job_registry(tmp_path,monkeypatch):
    monkeypatch.setenv('NANOPORE_JOB_DIR',str(tmp_path/'job-registry'))
