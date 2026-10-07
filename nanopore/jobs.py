"""Persistent local jobs independent of Streamlit/browser sessions."""
import json
import os
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


def job_root():
    configured=os.environ.get('NANOPORE_JOB_DIR')
    if configured:
        root=Path(configured).expanduser()
    elif Path('/workspace') in Path.cwd().resolve().parents or Path.cwd().resolve()==Path('/workspace'):
        root=Path('/workspace/.onboarding/nanopore-jobs')
    else:
        root=Path.home()/'nanopore-jobs'
    root.mkdir(parents=True,exist_ok=True)
    return root.resolve()


def atomic_json(path,data):
    path=Path(path)
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,default=str))
    temporary.replace(path)


def job_path(job_id):
    if len(job_id)!=32 or any(c not in '0123456789abcdef' for c in job_id):
        raise ValueError('Invalid job ID.')
    path=job_root()/job_id
    if not path.is_dir(): raise ValueError('Job not found on this computer.')
    return path


def process_token(pid):
    try:
        # Linux / WSL process start time prevents a reused PID appearing active.
        text=Path(f'/proc/{pid}/stat').read_text()
        fields=text[text.rfind(')')+2:].split()
        if fields[0]=='Z': return None
        boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        return boot+':'+fields[19]
    except (OSError,IndexError):
        try: os.kill(pid,0)
        except (OSError,ValueError): return None
        return 'present'


def status(job_id):
    path=job_path(job_id)
    data=json.loads((path/'status.json').read_text())
    if data['state']=='queued' and not data.get('pid') and (datetime.now(timezone.utc)-datetime.fromisoformat(data['created_at'])).total_seconds()>30:
        data.update(state='interrupted',message='Launch was interrupted before the worker started. Submit a new job.')
    if data['state'] in ('queued','running') and data.get('pid'):
        current=process_token(data['pid'])
        if not current or current!=data.get('process_token'):
            data['state']='interrupted'
            data['message']='Worker exited without completing. Inspect worker.log; partial outputs are not complete results.'
            # Read-only observation avoids races with worker completion.
    return data


def list_jobs(kind=None):
    jobs=[]
    for path in job_root().iterdir():
        if not path.is_dir() or not (path/'status.json').exists(): continue
        try:
            data=status(path.name)
            if kind is None or data.get('kind')==kind: jobs.append(data)
        except (ValueError,OSError,json.JSONDecodeError): continue
    return sorted(jobs,key=lambda x:x['created_at'],reverse=True)


def submit_job(request,uploads=None,table=None):
    """Stage inputs, then detach a worker. Local read files are referenced, not copied."""
    # This app targets Linux/macOS/WSL. Serialize competing submitters locally.
    import fcntl
    lock=job_root()/'.submit.lock'
    with lock.open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        if any(x['state'] in ('queued','running') for x in list_jobs()):
            raise ValueError('Another background analysis is active. You can browse its status; wait for it to finish before starting another heavy job.')
        job_id=uuid.uuid4().hex
        path=job_root()/job_id
        path.mkdir(mode=0o700)
        if uploads:
            import hashlib
            samples={}
            provenance={}
            for sample,files in uploads.items():
                if not files: continue
                directory=path/'inputs'/sample
                directory.mkdir(parents=True)
                samples[sample]=[]
                provenance[sample]=[]
                for i,file in enumerate(files):
                    name=Path(file.name.replace('\\','/')).name
                    target=directory/f'{i+1}_{name}'
                    with file.getbuffer() as view, target.open('wb') as output:
                        output.write(view)
                        digest=hashlib.sha256(view).hexdigest()
                    provenance[sample].append(dict(name=name,sha256=digest))
                    samples[sample].append(str(target))
            request={**request,'samples':samples,'input_provenance':provenance}
        if table is not None: save_frame(path/'input_table.json',table)
        atomic_json(path/'request.json',request)
        created=datetime.now(timezone.utc).isoformat()
        data=dict(id=job_id,kind=request['kind'],state='queued',message='Inputs saved; launching worker…',created_at=created)
        atomic_json(path/'status.json',data)
        try:
            with (path/'worker.log').open('ab') as log:
                process=subprocess.Popen([sys.executable,'-m','nanopore.worker',str(path)],cwd=str(Path(__file__).resolve().parents[1]),
                    stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True,close_fds=True)
            # Worker waits for this launch record before updating status.
            data.update(pid=process.pid,process_token=process_token(process.pid))
            atomic_json(path/'status.json',data)
            atomic_json(path/'launched.json',dict(pid=process.pid))
        except Exception as exc:
            atomic_json(path/'status.json',{**data,'state':'failed','message':str(exc)})
            raise
    return job_id


def save_frame(path,frame):
    frame.to_json(path,orient='table',index=False)


def load_frame(path):
    import pandas as pd
    return pd.read_json(path,orient='table')


def load_job(job_id):
    path=job_path(job_id)
    if status(job_id)['state']!='completed':
        raise ValueError('Results are available only after the worker completes successfully.')
    request=json.loads((path/'request.json').read_text())
    if request['kind']=='reads':
        result={p.stem:load_frame(p) for p in (path/'results').glob('*.json')}
        metadata=json.loads((path/'metadata.json').read_text())
        return result,metadata
    summary=load_frame(path/'summary.json')
    members=load_frame(path/'members.json')
    frame=load_frame(path/'input_table.json')
    dataset=(frame,request['reference'],request['annotations'],request['trim'],request['identity'],request['span'],request['context_identity'])
    return (summary,members,len(frame)),dataset
