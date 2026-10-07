"""Detached analysis process; never imports Streamlit or relies on its session."""
import json
import sys
import time
import traceback
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from .jobs import atomic_json, save_frame, load_frame
from .local import analyze_local
from .summary import summarize_table
from .workflow import runtime_info


def run_job(path):
    path=Path(path)
    for _ in range(100):
        if (path/'launched.json').exists(): break
        time.sleep(.05)
    else: raise RuntimeError('Worker launch record was not written.')
    state=json.loads((path/'status.json').read_text())
    request=json.loads((path/'request.json').read_text())
    last_update=0
    def update(message,force=False):
        nonlocal last_update
        now=time.monotonic()
        if not force and now-last_update<1: return
        state.update(message=message,updated_at=datetime.now(ZoneInfo('Europe/Stockholm')).isoformat())
        atomic_json(path/'status.json',state)
        last_update=now
    try:
        state['state']='running'
        update('Analysis running in a detached process. You may close the browser tab.',True)
        if request['kind']=='reads':
            samples={name:[Path(p) for p in files] for name,files in request['samples'].items()}
            result,report_metadata=analyze_local(samples,request['settings'],request['references'],request['output_root'],update)
            metadata=dict(settings=request['settings'],tools=request['tools'],references=request['references'],**report_metadata,
                background_job_id=state['id'],files={name:[dict(name=str(p),size_bytes=p.stat().st_size,modified_ns=p.stat().st_mtime_ns) for p in files] for name,files in samples.items()},
                workflow_runtime=runtime_info(),workflow_completed_at=datetime.now(ZoneInfo('Europe/Stockholm')).isoformat())
            if request.get('input_provenance'):
                metadata['files']=request['input_provenance']
            directory=path/'results';directory.mkdir()
            for name,frame in result.items(): save_frame(directory/f'{name}.json',frame)
            atomic_json(path/'metadata.json',metadata)
            atomic_json(Path(metadata['report_directory'])/'analysis_settings.json',metadata)
        elif request['kind']=='dal2':
            frame=load_frame(path/'input_table.json')
            summary,members=summarize_table(frame,request['reference'],request['annotations'],request['trim'],request['identity'],request['span'],
                lambda fraction:update(f'Analyzed {round(fraction*len(frame)):,} / {len(frame):,} Excel/CSV rows'),request['context_identity'])
            save_frame(path/'summary.json',summary);save_frame(path/'members.json',members)
            summary.to_csv(path/'dal2-retained-start-summary.csv',index=False)
            members.to_csv(path/'dal2-group-members.csv',index=False)
        else: raise ValueError('Unknown analysis job kind.')
        state['state']='completed'
        update('Analysis completed. Saved results are ready to load.',True)
    except BaseException as exc:
        state['state']='failed'
        update(f'Analysis failed: {type(exc).__name__}: {exc}. Partial outputs must not be treated as complete.',True)
        traceback.print_exc()
        return 1
    return 0


if __name__=='__main__':
    sys.exit(run_job(sys.argv[1]))
