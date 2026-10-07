import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from nanopore.io import references
from nanopore.dal2 import REFERENCE, ANNOTATIONS, FORWARD_PRIMER, REVERSE_PRIMER
from nanopore.alignment import load_consensus, build_rows, report_frame
from nanopore.viewer import viewer_html
from nanopore.summary import summarize_table
import plotly.express as px

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
context_identity=st.slider('Minimum downstream context identity for start grouping (%)',50,100,80)/100
st.caption('Grouping checks the next 50 reference bases from each annotated ATG (at least 40 covered bases), plus a 12-base near-start check (at least 10 covered bases and at least 85% identity). A different 5′ prefix is allowed and does not count as an earlier uORF. Weak downstream support and altered well-supported start codons remain separately flagged.')
frame=None
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
        if frame is None: raise ValueError('Upload a grouped consensus table first.')
        reference,annotations=alignment_inputs()
        progress=st.progress(0,text=f'Analyzing all {len(frame):,} imported rows…')
        def update_progress(fraction):
            progress.progress(fraction,text=f'Analyzed {round(fraction*len(frame)):,} / {len(frame):,} rows')
        summary,members=summarize_table(frame,reference,annotations,trim,identity,span,update_progress,context_identity)
        progress.empty()
        st.session_state['dal2_summary']=(summary,members,len(frame))
        st.session_state['dal2_dataset']=(frame.copy(),reference,annotations,trim,identity,span,context_identity)
        st.session_state.pop('dal2_view',None)
    except Exception as exc:
        st.error(f'Analysis could not finish: {exc}. No rows from this run have been silently skipped.')

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
