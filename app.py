import json
import hashlib
import streamlit as st
import pandas as pd
import plotly.express as px
from nanopore.io import dna, load_reads, references, fasta
from nanopore.pipeline import analyze
from nanopore.tools import dependencies
from nanopore.local import analyze_local, local_path
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from nanopore.workflow import workflow_excel, runtime_info

st.set_page_config(page_title='Nanopore 5′ RACE Amplicon Analyzer', page_icon='🧬', layout='wide')
st.title('Nanopore 5′ RACE Amplicon Analyzer')
st.page_link('pages/1_DAL2_Alignment_Viewer.py',label='Open DAL2 alignment viewer — import grouped CSV/Excel results',icon='🧬')
st.write('Explore amplicons or match known references in three independently sequenced PCR samples.')
st.warning('Read percentages describe PCR-derived sequencing-read abundance, not original RNA molecule abundance. Preliminary clusters and consensus sequences do not establish biological isoforms. PCR chimeras and sequencing artifacts require further validation.')
with st.expander('Before you begin · formats, primers and interpretation'):
    st.markdown('''**Input:** FASTQ, FASTQ.gz or FASTA; multiple files per sample are combined within that sample only. Reads from separate samples remain separate. **POD5 needs basecalling with Dorado before FASTQ analysis.**

**Primers:** enter each actual oligo in its **5′→3′ direction**. The anchor should occur upstream on the normalized strand; the reverse complement of the gene-specific primer should occur downstream. Both read orientations are searched. Mismatch allowance is a total edit budget per primer, including substitutions, insertions and deletions. Ambiguous IUPAC bases are supported in primers.

**Consensus:** VSEARCH generates an unpolished preliminary consensus. Single-primer reads are incomplete evidence. Similar-length sequences are never assumed identical. Cross-sample comparison uses only exact consensus DNA matches; sequencing errors may split the same biological sequence.

**References:** minimap2 uses map-ont and searches both strands. Identity is matching alignment bases divided by alignment block length (including gaps); coverage is aligned reference span divided by reference length. Separate alignments are not stitched. All qualifying multi-reference reads are marked ambiguous and excluded from unique reference totals. Minimap2 is heuristic; very short or highly divergent targets may not align.''')

with st.sidebar:
    st.header('Analysis settings')
    mode=st.radio('Analysis mode',['Amplicon discovery','Reference matching'],key='analysis_mode')
    min_quality=st.number_input('Minimum read Q score',0.0,60.0,7.0,.5,help='Mean read Q calculated from the mean base error probability.')
    min_length=st.number_input('Minimum read length (bp)',1,10000000,100)
    max_length=st.number_input('Maximum read length (bp)',1,10000000,20000)
    retain_fasta=st.checkbox('Include FASTA reads with unknown quality',value=True,help='FASTA has no Q scores. These reads pass the quality step only when this is checked.')
    if mode=='Amplicon discovery':
        edits=st.slider('Maximum edits per primer',0,8,2)
        cluster_identity=st.slider('Clustering sequence identity (%)',70.0,100.0,90.0,.5)/100
        keep_single=st.checkbox('Retain single-primer reads',value=False)
        full=st.checkbox('Keep primer-containing amplicons',value=False,help='Both primers: from anchor start through gene-primer site end. Single primer: retain the whole oriented read. Otherwise extract the interior or single-primer flank.')
    else:
        reference_identity=st.slider('Minimum alignment identity (%)',50.0,100.0,85.0,.5)/100
        reference_coverage=st.slider('Minimum reference coverage (%)',1.0,100.0,50.0,1.0)/100
        full_coverage=st.slider('Full-length reference coverage threshold (%)',50.0,100.0,95.0,1.0)/100
    with st.expander('Dependency checks',expanded=True):
        status=dependencies()
        for tool,value in status.items(): st.write(f'**{tool}:** {value}')
        st.caption('Install missing tools on your computer and restart Streamlit. See README for commands.')

st.subheader('1. Add sample files')
input_mode=st.radio('How to open reads', ['Upload small files', 'Read large files from this computer'],key='input_mode')
local_mode=input_mode=='Read large files from this computer'
if local_mode:
    st.info('No browser upload needed. Paste a full file path for each sample; one path per line. Files and samples are processed sequentially. Windows paths such as C:\\Users\\YourName\\Desktop\\reads.fastq.gz are accepted in Ubuntu/WSL. Leave unused samples empty.')
    st.caption('Keep .fastq.gz files compressed. VSEARCH still needs RAM for extracted amplicons from one sample; large or diverse samples may require substantial memory and time. Ubuntu/WSL may have a lower memory limit than your Windows installed RAM.')
    output_root=st.text_input('Folder to save complete reports',value=str(Path.home()/'nanopore-results'),help='Each run creates a separate folder. Reports remain on disk rather than being loaded into the browser.')
columns=st.columns(3)
uploads={}
local_samples={}
for i,col in enumerate(columns,1):
    with col:
        st.markdown(f'**Sample {i}**')
        if local_mode:
            paths=st.text_area('Full file paths (one per line)',key=f'paths_{i}',placeholder='C:\\Users\\YourName\\Desktop\\reads.fastq.gz')
            local_samples[f'Sample_{i}']=[local_path(line) for line in paths.splitlines() if line.strip()]
        else:
            uploads[f'Sample_{i}']=st.file_uploader('Sequencing files',type=['fastq','fq','gz','fasta','fa','fna'],accept_multiple_files=True,key=f'files_{i}')
st.caption('You can start with one sample; add all three for comparison. Files are processed in this local Streamlit process, never uploaded to an external analysis service. Larger datasets need more memory; prefilter or subsample externally if necessary.')

anchor=gene=''; ref_text=''
if mode=='Amplicon discovery':
    st.subheader('2. Enter primer sequences')
    a,b=st.columns(2)
    anchor=a.text_input('5′ RACE adapter / anchor primer (5′→3′)',placeholder='ACGT…')
    gene=b.text_input('Gene-specific primer (5′→3′)',placeholder='ACGT…')
else:
    st.subheader('2. Add reference sequences')
    ref_file=st.file_uploader('Reference FASTA',type=['fasta','fa','fna'],key='refs')
    ref_text=st.text_area('Or paste FASTA / a single DNA sequence',height=140)
    st.caption('If both are supplied, references are combined. Reference names must be unique.')

settings=dict(mode=mode,min_quality=min_quality,min_length=min_length,max_length=max_length,retain_fasta=retain_fasta)
if mode=='Amplicon discovery': settings.update(edits=edits,cluster_identity=cluster_identity,keep_single=keep_single,full=full,anchor=anchor,gene=gene)
else: settings.update(reference_identity=reference_identity,reference_coverage=reference_coverage,full_coverage=full_coverage)
# Results are explicitly a snapshot; do not silently relabel old results with changed controls.
if st.button('Analyze samples',type='primary'):
    try:
        if min_length>max_length: raise ValueError('Minimum read length must not exceed maximum read length.')
        if not any((local_samples if local_mode else uploads).values()): raise ValueError('Add sequencing files to at least one sample.')
        needed='vsearch' if mode=='Amplicon discovery' else 'minimap2'
        if status[needed].startswith(('Missing','Unable')): raise ValueError(f'{needed} is unavailable. Install it using the README instructions.')
        refs={}
        if mode=='Amplicon discovery':
            settings.update(anchor=dna(anchor),gene=dna(gene))
            if min(len(settings['anchor']),len(settings['gene']))<8: raise ValueError('Use primers at least 8 bases long.')
            if edits>=min(len(settings['anchor']),len(settings['gene']))/2: raise ValueError('Use an edit allowance below half the shortest primer length.')
        else:
            if full_coverage<reference_coverage: raise ValueError('Full-length coverage threshold must be at least the minimum reference coverage.')
            refs=references(ref_text)
            if ref_file:
                file_refs=references(ref_file.getvalue().decode('utf-8-sig'))
                if set(file_refs)&set(refs): raise ValueError('Reference names are duplicated between pasted and uploaded references.')
                refs.update(file_refs)
            if not refs: raise ValueError('Provide at least one reference sequence.')
        with st.status('Running analysis…',expanded=True) as task:
            if local_mode:
                if not output_root.strip(): raise ValueError('Choose a report output folder.')
                selected={name:paths for name,paths in local_samples.items() if paths}
                results,report_metadata=analyze_local(selected,settings,refs,local_path(output_root),lambda message: task.update(label=message))
                metadata=dict(settings=settings,tools=status,references=refs,**report_metadata,
                    files={name:[dict(name=str(path),size_bytes=path.stat().st_size,modified_ns=path.stat().st_mtime_ns) for path in paths] for name,paths in selected.items()})
                metadata.update(workflow_runtime=runtime_info(),workflow_completed_at=datetime.now(ZoneInfo('Europe/Stockholm')).isoformat())
                Path(metadata['report_directory'],'analysis_settings.json').write_text(json.dumps(metadata,indent=2))
            else:
                samples={name:load_reads(files,name) for name,files in uploads.items() if files}
                results=analyze(samples,settings,refs,task.write)
                metadata=dict(settings=settings,tools=status,files={name:[dict(name=f.name,sha256=hashlib.sha256(f.getvalue()).hexdigest()) for f in files] for name,files in uploads.items() if files},references=refs)
            if 'workflow_runtime' not in metadata:
                metadata.update(workflow_runtime=runtime_info(),workflow_completed_at=datetime.now(ZoneInfo('Europe/Stockholm')).isoformat())
            st.session_state['analysis']=(results,metadata)
            task.update(label='Analysis completed',state='complete',expanded=False)
    except Exception as exc:
        st.error(f'Analysis could not finish: {exc}')
        st.info('Check file format, primer orientation and thresholds. No partial results from this run are presented.')

if 'analysis' in st.session_state:
    result,metadata=st.session_state['analysis']
    st.subheader('Results')
    st.caption(f"Saved result snapshot · {metadata['settings']['mode']}. Changing inputs or settings requires clicking Analyze samples again.")
    st.dataframe(result['stats'],hide_index=True,use_container_width=True)
    if result['stats']['quality_filtered_reads'].sum()==0: st.warning('No reads passed filtering. Review length and Q thresholds; FASTA quality is unknown.')
    if metadata.get('local_files'):
        st.success('Full reports saved to: '+metadata['report_directory'])
        st.caption('Local-file mode: all reads are analyzed. Charts use weighted aggregate counts. Amplicon tables and plots show at most the first 500 clusters per sample; complete CSV/FASTA reports are on disk. Per-read reports are not loaded into the browser.')
        if st.button('Open report folder in Windows Explorer'):
            import subprocess
            import shutil
            if shutil.which('explorer.exe') and shutil.which('wslpath'):
                windows_path=subprocess.run(['wslpath','-w',metadata['report_directory']],capture_output=True,text=True,check=True).stdout.strip()
                subprocess.Popen(['explorer.exe',windows_path])
            else:
                st.info('Open the report folder using your file manager.')
    lengths=result['lengths']
    weighted={'y':'read_count','histfunc':'sum'} if 'read_count' in lengths else {}
    st.plotly_chart(px.histogram(lengths,**weighted,x='length',color='filter_status',facet_col='sample',nbins=60,title='Read-length distributions before and after filtering',labels={'length':'Read length (bp)'}),use_container_width=True)
    with st.expander('Quality statistics'):
        q=result.get('quality_bins',lengths).dropna(subset=['mean_q'])
        if len(q): st.plotly_chart(px.histogram(q,**({'y':'read_count','histfunc':'sum'} if 'read_count' in q else {}),x='mean_q',color='sample',nbins=40,title='Read Q scores (FASTQ only)'),use_container_width=True)
        else: st.info('No quality scores available; these are FASTA reads.')
    if metadata['settings']['mode']=='Amplicon discovery':
        table=result['amplicons']
        st.caption('percent_quality_filtered = cluster reads / all quality-filtered reads. percent_assigned = cluster reads / reads assigned to any amplicon. Single-primer reads count as assigned only when enabled. Clustering uses sequence identity, never length alone.')
        if table.empty: st.warning('No amplicons found. Check primer sequences, orientation, edits and filtering.')
        else:
            st.dataframe(table,hide_index=True,use_container_width=True)
            st.plotly_chart(px.scatter(table,x='length',y='supporting_reads',color='sample',hover_data=['amplicon','percent_quality_filtered'],title='Amplicon length versus read support',labels={'length':'Preliminary consensus length (bp)'}),use_container_width=True)
            sequences=sorted(set(table['consensus_dna']))
            labels={seq:f'Sequence {i+1}' for i,seq in enumerate(sequences)}
            compare=table.assign(sequence_group=table['consensus_dna'].map(labels)).groupby(['sample','sequence_group'],as_index=False)['percent_quality_filtered'].sum()
            st.plotly_chart(px.bar(compare,x='sequence_group',y='percent_quality_filtered',color='sample',barmode='group',title='Across-sample abundance · exact consensus DNA groups'),use_container_width=True)
            st.caption('Groups share exactly the same consensus DNA. Different sequences of the same length remain separate; this is not isoform classification.')
            st.download_button('Download preliminary amplicon FASTA'+(' (preview)' if metadata.get('local_files') else ''),fasta((f"{r.amplicon}|reads={r.supporting_reads}|preliminary",r.consensus_dna) for r in table.itertuples()),'preliminary_amplicons.fasta','text/plain')
    else:
        table=result['matches']
        st.caption('Unique assignments only. percent_quality_filtered uses all quality-filtered reads; percent_reference_assigned uses only uniquely reference-assigned reads. Ambiguous reads are excluded from both numerators and the assigned denominator. A full-length match meets the configured reference coverage threshold; it does not prove a full transcript or intact primer sites.')
        st.dataframe(table,hide_index=True,use_container_width=True)
        st.plotly_chart(px.bar(table,x='reference',y='percent_quality_filtered',color='sample',barmode='group',hover_data=['full_length_reads','partial_reads','ambiguous_candidate_reads'],title='Reference matches as % of all quality-filtered reads'),use_container_width=True)
        with st.expander('Per-read alignments and ambiguous candidates'):
            st.dataframe(result['match_details'],hide_index=True,use_container_width=True)
        st.download_button('Download reference FASTA',fasta(metadata['references'].items()),'references.fasta','text/plain')
    st.subheader('Download reports')
    st.download_button('Download workflow Excel (no sample results)',workflow_excel(metadata,settings_source='Completed read-analysis snapshot'),'workflow-and-settings.xlsx','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',key='workflow_read_download')
    for name,frame in result.items():
        if not frame.empty:
            st.download_button(f'Download {name.replace("_"," ")} CSV'+(' (preview/aggregated)' if metadata.get('local_files') else ''),frame.to_csv(index=False),f'{name}.csv','text/csv',key=f'dl_{name}')
    st.download_button('Download analysis settings and file provenance',json.dumps(metadata,indent=2),'analysis_settings.json','application/json')
