import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from nanopore.io import references
from nanopore.dal2 import REFERENCE, ANNOTATIONS, FORWARD_PRIMER, REVERSE_PRIMER
from nanopore.alignment import load_consensus, build_rows, report_frame
from nanopore.viewer import viewer_html

st.set_page_config(page_title='DAL2 5′ RACE alignment viewer',page_icon='🧬',layout='wide')
st.title('DAL2 5′ RACE alignment viewer')
st.write('Import grouped consensus sequences, align their shared DAL2 region, and arrange rows to compare observed 5′ ends.')
st.warning('Preliminary consensus alignment. Read abundance is not RNA abundance. Observed RACE boundaries are not proven transcription start sites; ATG retention is not evidence of translation or an intact complete ORF.')
with st.expander('Your experimental design and alignment interpretation',expanded=True):
    st.write('Default reference: 589 bases, starting at uORF1. Main ATG: 100; internal ATGs: 122 and 170. Reverse-primer binding site: 564–589. All coordinates are 1-based and inclusive.')
    st.code('PCR forward primer: '+FORWARD_PRIMER+'\nReverse primer: '+REVERSE_PRIMER+'\nResidual TSO after PCR forward primer: GCAGAGTACATGGG')
    st.write('The residual TSO can remain after the original app trims the PCR primer. Candidate adapter matches near the 5′ end can be removed for alignment, with the removed sequence reported. RNA rG residues are represented as G in DNA. The TSO/biological G boundary remains uncertain. Unrecognized or truncated adapter sequence is not automatically removed.')
    st.write('Local alignment allows mismatches and indels and tests both orientations. Upstream sequence outside the reference is preserved. Unaligned prefixes are displayed separately, not treated as matching each other. A downstream alignment with an unaligned prefix has an uncertain 5′ endpoint. Feature states describe ATG sites, not complete uORFs.')

upload=st.file_uploader('Grouped consensus CSV or Excel (.xlsx)',type=['csv','xlsx'])
reference_text=st.text_area('Reference DNA or single-record FASTA',value=REFERENCE,height=130)
with st.expander('Reference annotations — update these if you change the reference'):
    edited=st.data_editor(pd.DataFrame(ANNOTATIONS),num_rows='dynamic',hide_index=True,key='annotations',column_config={'kind':st.column_config.SelectboxColumn(options=['start','primer','region'])})
trim=st.checkbox('Remove candidate TSO-derived 5′ sequence before alignment',value=True)
a,b=st.columns(2)
identity=a.slider('Minimum aligned sequence identity (%)',50,100,70)/100
span=b.number_input('Minimum aligned reference span (bp)',20,20000,80)
frame=None
if upload:
    try:
        frame=load_consensus(upload.name,upload.getvalue())
        st.caption(f'{len(frame):,} consensus sequences available. Choose up to 100 rows per interactive comparison; all selected rows are reported, including uncertain alignments.')
        labels={i:f"{r.get('sample','')} · {r['amplicon']} · reads={r.get('supporting_reads','unknown')}" for i,r in enumerate(frame.to_dict('records'))}
        selected=st.multiselect('Amplicons to compare',options=list(labels),default=list(labels)[:min(30,len(labels))],format_func=lambda i:labels[i],max_selections=100)
    except Exception as exc:
        st.error(str(exc))
if st.button('Align selected amplicons',type='primary'):
    try:
        if frame is None: raise ValueError('Upload a grouped consensus table first.')
        if not selected: raise ValueError('Select at least one amplicon.')
        parsed=references(reference_text)
        if len(parsed)!=1: raise ValueError('Provide exactly one reference for this annotated viewer.')
        reference=next(iter(parsed.values()))
        if len(reference)>10000: raise ValueError('Use a reference of at most 10,000 bases for this viewer.')
        annotations=edited.dropna(how='all').to_dict('records')
        for annotation in annotations:
            annotation['start']=int(annotation['start']);annotation['end']=int(annotation['end'])
        annotations.sort(key=lambda x:x['start'])
        progress=st.progress(0,text='Aligning preliminary consensus sequences…')
        rows=build_rows(frame.iloc[selected],reference,annotations,trim,identity,span,progress.progress)
        progress.empty()
        st.session_state['dal2_view']=(reference,annotations,rows)
    except Exception as exc:
        st.error(f'Alignment could not finish: {exc}')

if 'dal2_view' in st.session_state:
    reference,annotations,rows=st.session_state['dal2_view']
    st.subheader('Interactive alignment — drag rows to reorder')
    st.caption('This is a saved alignment snapshot; changing inputs requires Align selected amplicons again. Row order is managed in the viewer: use its SVG and row-order exports to preserve your arrangement.')
    html=viewer_html(reference,annotations,rows)
    components.html(html,height=980,scrolling=True)
    st.caption('Reference-projection FASTA includes aligned reference columns only: insertions and unaligned extensions are omitted. Use the full CSV/JSON for those sequences and the self-contained HTML for the complete interactive display.')
    report=report_frame(rows)
    st.dataframe(report.drop(columns=['sequence','upstream_sequence','downstream_sequence','adapter_removed'],errors='ignore'),hide_index=True)
    st.download_button('Download alignment and feature CSV',report.to_csv(index=False),'dal2-alignment-report.csv','text/csv')
    st.download_button('Download complete alignment JSON',json.dumps(dict(reference=reference,annotations=annotations,rows=rows),indent=2,default=str),'dal2-alignment.json','application/json')
    st.download_button('Download standalone interactive viewer',html,'dal2-interactive-alignment.html','text/html')
