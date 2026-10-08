import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from nanopore.io import references
from nanopore.dal2 import REFERENCE, ANNOTATIONS, FORWARD_PRIMER, REVERSE_PRIMER
from nanopore.alignment import load_consensus, build_rows, report_frame
from nanopore.viewer import viewer_html
from nanopore.workflow import workflow_excel, read_settings_metadata
import plotly.express as px
from nanopore.jobs import submit_job
from nanopore.job_ui import job_panel
from nanopore.local import local_path

st.set_page_config(page_title='DAL2 5′ RACE alignment viewer',page_icon='🧬',layout='wide')
st.title('DAL2 5′ RACE alignment viewer')
st.write('Import grouped consensus sequences, align their shared DAL2 region, and arrange rows to compare observed 5′ ends.')
st.warning('Preliminary consensus alignment. Read abundance is not RNA abundance. Observed RACE boundaries are not proven transcription start sites; ATG retention is not evidence of translation or an intact complete ORF.')
with st.expander('Your experimental design and alignment interpretation',expanded=True):
    st.write('Default reference: 589 bases, starting at uORF1. Main ATG: 100; internal ATGs: 122 and 170. Reverse-primer binding site: 564–589. All coordinates are 1-based and inclusive.')
    st.code('PCR forward primer: '+FORWARD_PRIMER+'\nReverse primer: '+REVERSE_PRIMER+'\nResidual TSO after PCR forward primer: GCAGAGTACATGGG')
    st.write('The residual TSO can remain after the original app trims the PCR primer. Candidate adapter matches near the 5′ end can be removed for alignment, with the removed sequence reported. RNA rG residues are represented as G in DNA. The TSO/biological G boundary remains uncertain. Unrecognized or truncated adapter sequence is not automatically removed.')
    st.write('Local alignment allows mismatches and indels and tests both orientations. Upstream sequence outside the reference is preserved. Unaligned prefixes are displayed separately, not treated as matching each other. A downstream alignment with an unaligned prefix has an uncertain 5′ endpoint. Feature states describe ATG sites, not complete uORFs.')

input_mode=st.radio('Grouped consensus input', ['Upload a file', 'Use a file on this computer'], horizontal=True)
upload=None
frame=None
if input_mode=='Upload a file':
    upload=st.file_uploader('Grouped consensus CSV or Excel (.xlsx)',type=['csv','xlsx'])
else:
    st.caption('Bypasses the browser upload size limit. In Windows Explorer, right-click your CSV or Excel file and choose Copy as path, then paste it below. Windows paths are converted automatically for Ubuntu/WSL. The file must be on the computer running this app.')
    path_text=st.text_input('Full path to grouped consensus CSV or Excel',placeholder=r'C:\Users\Siguradottir\Downloads\amplicons.csv')
    if st.button('Load file from this computer'):
        st.session_state.pop('dal2_local_input',None)
        try:
            if not path_text.strip(): raise ValueError('Paste a file path first.')
            path=local_path(path_text)
            if not path.is_file(): raise ValueError('File not found. Use Copy as path on the file itself, not its folder.')
            with st.spinner('Reading and validating every row…'):
                loaded=load_consensus(path)
            st.session_state['dal2_local_input']=(str(path),loaded)
        except Exception as exc:
            st.error(f'Cannot load this file: {exc}')
    saved=st.session_state.get('dal2_local_input')
    if saved and saved[0]==str(local_path(path_text)):
        frame=saved[1]
        st.info(f'{len(frame):,} rows imported from {saved[0]}. Analysis will process EVERY row. Click Load again if you edit the source file.')
reference_text=st.text_area('Reference DNA or single-record FASTA',value=REFERENCE,height=130)
with st.expander('Reference annotations — update these if you change the reference'):
    edited=st.data_editor(pd.DataFrame(ANNOTATIONS),num_rows='dynamic',hide_index=True,key='annotations',column_config={'kind':st.column_config.SelectboxColumn(options=['start','primer','region'])})
trim=st.checkbox('Remove candidate TSO-derived 5′ sequence before alignment',value=True)
a,b=st.columns(2)
identity=a.slider('Minimum aligned sequence identity (%)',50,100,70)/100
span=b.number_input('Minimum aligned reference span (bp)',20,20000,80)
context_identity=st.slider('Minimum downstream context identity for start grouping (%)',50,100,80)/100
st.caption('Grouping checks the next 50 reference bases from each annotated ATG (at least 40 covered bases), plus a 12-base near-start check (at least 10 covered bases and at least 85% identity). A different 5′ prefix is allowed and does not count as an earlier uORF. Weak downstream support and altered well-supported start codons remain separately flagged.')
if upload:
    try:
        frame=load_consensus(upload.name,upload.getvalue())
        st.info(f'{len(frame):,} rows imported. Analysis will process EVERY row; the interactive display is a separate selection made afterward.')
    except Exception as exc:
        st.error(str(exc))

def alignment_inputs():
    parsed=references(reference_text)
    if len(parsed)!=1: raise ValueError('Provide exactly one reference for this annotated viewer.')
    reference=next(iter(parsed.values()))
    if len(reference)>10000: raise ValueError('Use a reference of at most 10,000 bases for this viewer.')
    annotations=edited.dropna(how='all').to_dict('records')
    for annotation in annotations:
        annotation['start']=int(annotation['start']);annotation['end']=int(annotation['end'])
    annotations.sort(key=lambda x:x['start'])
    return reference,annotations

if st.button('Analyze ALL Excel/CSV rows',type='primary'):
    try:
        if frame is None: raise ValueError('Upload a grouped consensus table or load a file from this computer first.')
        reference,annotations=alignment_inputs()
        request=dict(kind='dal2',reference=reference,annotations=annotations,trim=trim,identity=identity,span=span,context_identity=context_identity)
        with st.spinner('Saving the complete table and launching background analysis…'):
            job_id=submit_job(request,table=frame)
        st.session_state['active_job_dal2']=job_id
        st.session_state['job_select_dal2']=job_id
        st.success('Background table analysis started. You may close the browser tab and reconnect using the job list.')
    except Exception as exc:
        st.error(f'Analysis could not finish: {exc}. No rows from this run have been silently skipped.')

job_panel('dal2')

if 'dal2_summary' in st.session_state:
    summary,members,total=st.session_state['dal2_summary']
    st.subheader(f'Full-table retained-start summary · {total:,} imported amplicons')
    st.success(f'{total:,} imported rows · {len(members):,} analyzed rows · {int(summary.amplicon_clusters.sum()):,} rows accounted for in groups')
    st.caption('Saved analysis snapshot: rerun Analyze ALL after changing inputs or settings. Different/unaligned 5′ sequences do not block a downstream uORF group. Groups represent the first intact annotated ATG with supported downstream sequence, not exact biological starts or complete ORFs. Earlier reference uORFs outside the alignment are not inferred from the unaligned prefix.')
    st.dataframe(summary,hide_index=True)
    st.caption('percent_of_imported_supporting_reads uses all supporting reads in this imported table per sample. The other percentages retain the original quality-filtered/assigned denominators. A partial source file still gives only a subset summary; use the complete amplicons.csv. Missing counts remain unknown. Read counts assume disjoint clusters.')
    if summary.supporting_reads.notna().all():
        st.plotly_chart(px.bar(summary,x='retained_start_group',y='supporting_reads',color='sample',barmode='group',title='Read support by retained start · uORF order'),use_container_width=True)
    for sample in summary['sample'].unique():
        for group in summary.loc[summary['sample']==sample,'retained_start_group']:
            subset=members[(members['sample']==sample)&(members['retained_start_group']==group)]
            with st.expander(f'{sample} · {group} · {len(subset):,} amplicon clusters'):
                st.dataframe(subset.drop(columns=['sequence','upstream_sequence','downstream_sequence','adapter_removed'],errors='ignore'),hide_index=True)
                st.caption('All member sequences and different 5′ prefixes are preserved in the full downloadable CSV.')
    st.download_button('Download retained-start summary CSV',summary.to_csv(index=False),'dal2-retained-start-summary.csv','text/csv')
    st.download_button('Download ALL analyzed rows and sequences CSV',members.to_csv(index=False),'dal2-group-members.csv','text/csv')

    if 'dal2_dataset' in st.session_state:
        dataset,reference,annotations,saved_trim,saved_identity,saved_span,saved_context=st.session_state['dal2_dataset']
        with st.expander('Download shared workflow and settings Excel',expanded=True):
            st.caption('One workflow report, without sample names or counts. It uses this completed alignment snapshot. To include the earlier FASTQ settings after a restart, upload its saved analysis_settings.json. These settings are not independently matched to the imported consensus file.')
            settings_file=st.file_uploader('Earlier FASTQ analysis settings JSON (optional)',type=['json'],key='workflow_settings_json')
            read_metadata=None
            source='Not supplied: earlier FASTQ settings are not recorded'
            valid_settings=True
            if settings_file:
                try:
                    read_metadata=read_settings_metadata(settings_file.getvalue())
                    source='User-supplied analysis_settings.json; association with consensus input not independently verified'
                except Exception as exc:
                    st.error(f'Cannot use settings file: {exc}')
                    valid_settings=False
            elif 'analysis' in st.session_state:
                read_metadata=st.session_state['analysis'][1]
                source='Completed read-analysis snapshot in this session; association with consensus input not independently verified'
            if valid_settings:
                alignment_settings=dict(reference=reference,annotations=annotations,trim=saved_trim,identity=saved_identity,span=saved_span,context_identity=saved_context)
                st.download_button('Download workflow Excel (no sample results)',workflow_excel(read_metadata,alignment_settings,source),'workflow-and-settings.xlsx','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',key='workflow_alignment_download')
        st.subheader('Choose analyzed rows for the interactive image')
        group=st.selectbox('Display a start group',options=['All groups']+summary.retained_start_group.drop_duplicates().tolist())
        eligible=members if group=='All groups' else members[members.retained_start_group==group]
        group_rank={name:i for i,name in enumerate(summary.retained_start_group.drop_duplicates())}
        positions=[int(row_id.split('_')[-1])-1 for row_id in eligible.id]
        positions.sort(key=lambda i:group_rank[members.iloc[i]['retained_start_group']])
        labels={i:f"{members.iloc[i]['sample']} · {members.iloc[i]['name']} · {members.iloc[i]['retained_start_group']}" for i in positions}
        selected=st.multiselect('Rows to display (analysis already includes ALL rows)',options=positions,default=positions[:30],format_func=lambda i:labels[i],max_selections=100,key='display_rows_'+group)
        st.caption('The 100-row limit applies only to the interactive image, never the full analysis or CSV exports. Rebuilding the image uses the saved analysis settings.')
        if st.button('Build interactive image of selected analyzed rows'):
            if not selected:
                st.error('Choose at least one analyzed row to display.')
            else:
                progress=st.progress(0,text='Building interactive image…')
                rows=build_rows(dataset.iloc[selected],reference,annotations,saved_trim,saved_identity,saved_span,progress.progress)
                for row,idx in zip(rows,selected):
                    row['metadata']['retained_start_group']=members.iloc[idx]['retained_start_group']
                progress.empty()
                st.session_state['dal2_view']=(reference,annotations,rows)

if 'dal2_view' in st.session_state:
    reference,annotations,rows=st.session_state['dal2_view']
    st.subheader('Interactive alignment — drag rows to reorder')
    st.caption('Image selection only. The full-table summary and exports include every analyzed row. Row order is managed in the viewer; use its SVG and order exports to preserve your arrangement.')
    html=viewer_html(reference,annotations,rows)
    components.html(html,height=980,scrolling=True)
    st.caption('Reference-projection FASTA omits insertions and unaligned extensions. Use CSV/JSON for full sequences and HTML for the complete interactive display.')
    report=report_frame(rows)
    st.dataframe(report.drop(columns=['sequence','upstream_sequence','downstream_sequence','adapter_removed'],errors='ignore'),hide_index=True)
    st.download_button('Download displayed-row alignment CSV',report.to_csv(index=False),'dal2-displayed-alignment.csv','text/csv')
    st.download_button('Download displayed-row alignment JSON',json.dumps(dict(reference=reference,annotations=annotations,rows=rows),indent=2,default=str),'dal2-displayed-alignment.json','application/json')
    st.download_button('Download standalone interactive viewer',html,'dal2-interactive-alignment.html','text/html')
