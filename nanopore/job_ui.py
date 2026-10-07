"""Small status panel; disk job registry survives browser/Streamlit sessions."""
import streamlit as st
from datetime import datetime
from zoneinfo import ZoneInfo
from .jobs import list_jobs, status, load_job, job_path


def job_panel(kind):
    st.subheader('Background jobs')
    st.caption('Jobs continue when you close the browser or change pages. Keep the computer awake and Ubuntu/WSL running. Power-off, reboot or WSL shutdown interrupts workers; interrupted jobs must be restarted. Reopen this page to load saved results.')
    if st.button('Refresh job list',key=f'job_refresh_{kind}'): st.rerun()
    jobs=list_jobs(kind)
    if not jobs:
        st.info('No saved background jobs for this page yet.')
        return
    ids=[j['id'] for j in jobs]
    preferred=st.session_state.get(f'active_job_{kind}')
    def label(job_id):
        job=next(j for j in jobs if j['id']==job_id)
        time=datetime.fromisoformat(job['created_at']).astimezone(ZoneInfo('Europe/Stockholm')).strftime('%Y-%m-%d %H:%M')
        return f"{time} Stockholm · {job['state']} · {job_id[:8]}"
    selected=st.selectbox('Select background job',ids,index=ids.index(preferred) if preferred in ids else 0,
        format_func=label,key=f'job_select_{kind}')
    st.code('Job ID: '+selected+'\nSaved job folder: '+str(job_path(selected)))

    @st.fragment(run_every='3s')
    def live_status():
        data=status(selected)
        if data['state']=='completed': st.success(data['message'])
        elif data['state'] in ('failed','interrupted'): st.error(data['message'])
        else: st.info(data['message'])
        if data['state']=='completed' and st.button('Load completed results',key=f'job_load_{kind}_{selected}'):
            loaded=load_job(selected)
            if kind=='reads': st.session_state['analysis']=loaded
            else:
                st.session_state['dal2_summary'],st.session_state['dal2_dataset']=loaded
                st.session_state.pop('dal2_view',None)
            st.rerun(scope='app')
    live_status()
