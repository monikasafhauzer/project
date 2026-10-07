import json
import os
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd
import pytest
from nanopore.jobs import submit_job,status,load_job,list_jobs,job_path
from nanopore.dal2 import REFERENCE,ANNOTATIONS


def wait_for_job(job_id):
    deadline=time.monotonic()+20
    while time.monotonic()<deadline:
        state=status(job_id)
        if state['state'] not in ('queued','running'): return state
        time.sleep(.1)
    raise AssertionError('Background job did not finish within test timeout.')


def request():
    return dict(kind='dal2',reference=REFERENCE,annotations=ANNOTATIONS,trim=True,identity=.7,span=80,context_identity=.8)


def test_detached_worker_survives_submitter_exit(tmp_path):
    table=pd.DataFrame([dict(sample='S',amplicon='uorf2',consensus_dna='ACCTGACCTGACCTGA'+REFERENCE[13:563],supporting_reads=10)])
    source=tmp_path/'input.json';table.to_json(source,orient='table',index=False)
    config=tmp_path/'config.json';config.write_text(json.dumps(request()))
    # Simulates the browser/server-side submission returning and its process exiting.
    program='from nanopore.jobs import submit_job,load_frame;import json,sys;print(submit_job(json.load(open(sys.argv[1])),table=load_frame(sys.argv[2])),flush=True)'
    parent=subprocess.run([sys.executable,'-c',program,str(config),str(source)],capture_output=True,text=True,check=True)
    job_id=parent.stdout.strip()
    state=wait_for_job(job_id)
    assert state['state']=='completed',(state,(job_path(job_id)/'worker.log').read_text())
    summary,dataset=load_job(job_id)
    assert summary[0].supporting_reads.sum()==10
    assert summary[1].retained_start_group.iloc[0]=='uORF2 ATG retained'
    assert len(dataset[0])==1
    assert list_jobs('dal2')[0]['id']==job_id


def test_real_fastq_background_results(tmp_path):
    from test_local import write_fastq,settings
    from nanopore.io import reverse_complement
    seq=settings()['anchor']+REFERENCE[:563]+reverse_complement(settings()['gene'])
    path=tmp_path/'reads.fastq.gz';write_fastq(path,[seq,seq])
    config=dict(kind='reads',samples={'Sample_1':[str(path)]},settings=settings(),references={},output_root=str(tmp_path/'reports'),tools={})
    job_id=submit_job(config)
    state=wait_for_job(job_id)
    assert state['state']=='completed',(state,(job_path(job_id)/'worker.log').read_text())
    result,metadata=load_job(job_id)
    assert result['stats'].input_reads.iloc[0]==2
    assert result['amplicons'].supporting_reads.iloc[0]==2
    assert Path(metadata['report_directory'],'analysis_settings.json').exists()
    assert metadata['background_job_id']==job_id


def test_failure_is_saved_and_results_are_not_loadable(tmp_path):
    config=dict(kind='reads',samples={'S':[str(tmp_path/'missing.fastq.gz')]},settings={},references={},output_root=str(tmp_path/'reports'),tools={})
    job_id=submit_job(config)
    assert wait_for_job(job_id)['state']=='failed'
    with pytest.raises(ValueError,match='only after'): load_job(job_id)
    assert 'File not found' in status(job_id)['message']


def test_interrupted_status_and_invalid_ids(tmp_path):
    from nanopore.jobs import atomic_json,job_root
    root=job_root();path=root/('a'*32);path.mkdir()
    atomic_json(path/'status.json',dict(id=path.name,kind='reads',state='running',created_at='2026-10-07T00:00:00+00:00',pid=999999999,process_token='invalid',message='running'))
    assert status(path.name)['state']=='interrupted'
    with pytest.raises(ValueError): job_path('../file')


def test_second_active_job_is_rejected():
    from nanopore.jobs import atomic_json,job_root,process_token
    from datetime import datetime,timezone
    path=job_root()/('b'*32);path.mkdir()
    atomic_json(path/'status.json',dict(id=path.name,kind='reads',state='running',created_at=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),process_token=process_token(os.getpid()),message='running'))
    with pytest.raises(ValueError,match='Another background analysis'):
        submit_job(request(),table=pd.DataFrame())


def test_upload_staged_before_detaching(tmp_path):
    import io
    from test_local import settings
    class Upload(io.BytesIO):
        name='reads.fasta'
    seq=settings()['anchor']+REFERENCE[:563]+__import__('nanopore.io',fromlist=['reverse_complement']).reverse_complement(settings()['gene'])
    upload=Upload(('>r\n'+seq+'\n').encode())
    config=dict(kind='reads',samples={},settings=settings(),references={},output_root=str(tmp_path/'reports'),tools={})
    job_id=submit_job(config,uploads={'Sample_1':[upload]})
    upload.close()
    assert wait_for_job(job_id)['state']=='completed'
    result,_=load_job(job_id)
    assert result['stats'].input_reads.iloc[0]==1
