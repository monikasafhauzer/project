import pandas as pd
import pytest
from nanopore.dal2 import REFERENCE,ANNOTATIONS
from nanopore.summary import summarize_table,UNCERTAIN


def table():
    return pd.DataFrame([
        dict(sample='S1',amplicon='upstream_a',consensus_dna='ACCTGACCTGACCTGA'+REFERENCE[:563],supporting_reads=10,percent_quality_filtered=1,percent_assigned=2),
        dict(sample='S1',amplicon='upstream_b',consensus_dna='TCTACCTACCTACCTA'+REFERENCE[:563],supporting_reads=20,percent_quality_filtered=2,percent_assigned=4),
        dict(sample='S1',amplicon='uorf2',consensus_dna=REFERENCE[13:563],supporting_reads=30,percent_quality_filtered=3,percent_assigned=6),
        dict(sample='S1',amplicon='internal',consensus_dna=REFERENCE[169:563],supporting_reads=40,percent_quality_filtered=4,percent_assigned=8),
        dict(sample='S2',amplicon='internal',consensus_dna=REFERENCE[169:563],supporting_reads=7,percent_quality_filtered=10,percent_assigned=100)])


def test_full_table_grouping_and_sample_denominators():
    summary,members=summarize_table(table(),REFERENCE,ANNOTATIONS)
    s=summary[summary['sample']=='S1'].set_index('retained_start_group')
    assert s.loc['uORF1 ATG retained','amplicon_clusters']==2
    assert s.loc['uORF1 ATG retained','supporting_reads']==30
    assert s.loc['uORF1 ATG retained','percent_of_imported_supporting_reads']==30
    assert s.loc['uORF1 ATG retained','percent_quality_filtered_from_input']==3
    assert s.loc['uORF1 ATG retained','percent_assigned_from_input']==6
    assert s.loc['uORF1 ATG retained','clusters_with_unaligned_5prime_extension']==2
    assert s.supporting_reads.sum()==100
    assert s.percent_of_imported_supporting_reads.sum()==100
    assert len(members)==5 and 'bases' not in members
    assert summary[summary['sample']=='S2'].percent_of_imported_supporting_reads.iloc[0]==100


def test_nonreference_prefix_does_not_block_supported_start():
    frame=pd.DataFrame([
        dict(amplicon='divergent_prefix',consensus_dna='ACCTACCGACTACCTA'+REFERENCE[169:563],supporting_reads=10),
        dict(amplicon='mutated_uorf1',consensus_dna=REFERENCE[:1]+'C'+REFERENCE[2:563],supporting_reads=20),
        dict(amplicon='unmatched',consensus_dna='ACGT'*70,supporting_reads=30)])
    summary,members=summarize_table(frame,REFERENCE,ANNOTATIONS,min_identity=.9)
    assert members.retained_start_group.iloc[0]=='Internal ATG B retained'
    assert members.retained_start_group.iloc[2]==UNCERTAIN
    assert summary.supporting_reads.sum()==60


def test_missing_counts_and_duplicate_rows():
    frame=table();frame.loc[0,'supporting_reads']=None
    summary,_=summarize_table(frame,REFERENCE,ANNOTATIONS)
    assert summary[summary['sample']=='S1'].percent_of_imported_supporting_reads.isna().all()
    with pytest.raises(ValueError,match='Duplicate'):
        summarize_table(pd.concat([table(),table().iloc[:1]]),REFERENCE,ANNOTATIONS)
    for value in (-1,1.5,'abc'):
        bad=table();bad['supporting_reads']=bad.supporting_reads.astype(object);bad.loc[0,'supporting_reads']=value
        with pytest.raises(ValueError,match='supporting_reads'):
            summarize_table(bad,REFERENCE,ANNOTATIONS)


def test_summary_page_rendering():
    from pathlib import Path
    from streamlit.testing.v1 import AppTest
    page=Path(__file__).resolve().parents[1]/'pages'/'1_DAL2_Alignment_Viewer.py'
    summary,members=summarize_table(table(),REFERENCE,ANNOTATIONS)
    app=AppTest.from_file(str(page))
    app.session_state['dal2_summary']=(summary,members,5)
    app.run(timeout=30)
    assert not app.exception
    assert any('Full-table' in x.value for x in app.subheader)
    assert len(app.expander)>=6


def test_all_annotated_starts_with_different_upstream_sequences():
    from test_analysis import random_dna
    starts=[a for a in ANNOTATIONS if a['kind']=='start']
    records=[]
    for i,start in enumerate(starts):
        # An unmatched prefix can contain ATGs without being an earlier DAL2 uORF.
        records.append(dict(sample='S',amplicon=f'variant_{i}',consensus_dna=random_dna(60,100+i)+REFERENCE[start['start']-1:563],supporting_reads=i+1))
    summary,members=summarize_table(pd.DataFrame(records),REFERENCE,ANNOTATIONS)
    assert members.retained_start_group.tolist()==[a['label']+' retained' for a in starts]
    assert summary.retained_start_group.tolist()==[a['label']+' retained' for a in starts]
    assert summary.supporting_reads.sum()==sum(range(1,len(starts)+1))
    assert all('unaligned' in x for x in members.upstream_interpretation)


def test_full_analysis_beyond_old_row_limit(monkeypatch):
    import nanopore.summary as module
    from nanopore.alignment import align_consensus,load_consensus
    result=align_consensus(REFERENCE[13:563],REFERENCE,ANNOTATIONS)
    calls=[]
    def mock_alignment(*args):
        calls.append(1)
        return dict(result)
    monkeypatch.setattr(module,'align_consensus',mock_alignment)
    n=5001
    frame=pd.DataFrame([dict(amplicon=f'a{i}',consensus_dna=REFERENCE[13:563],supporting_reads=1) for i in range(n)])
    loaded=load_consensus('all.csv',frame.to_csv(index=False).encode())
    summary,members=summarize_table(loaded,REFERENCE,ANNOTATIONS)
    assert len(calls)==n and len(members)==n
    assert summary.amplicon_clusters.sum()==n and summary.supporting_reads.sum()==n


def test_excel_import_includes_all_data_worksheets():
    import io
    from nanopore.alignment import load_consensus
    buffer=io.BytesIO()
    with pd.ExcelWriter(buffer,engine='openpyxl') as writer:
        table().iloc[:2].to_excel(writer,sheet_name='First',index=False)
        table().iloc[2:].to_excel(writer,sheet_name='Second',index=False)
    frame=load_consensus('all.xlsx',buffer.getvalue())
    assert len(frame)==5 and set(frame.source_sheet)=={'First','Second'}


def test_interactive_selection_follows_complete_analysis():
    from pathlib import Path
    from streamlit.testing.v1 import AppTest
    page=Path(__file__).resolve().parents[1]/'pages'/'1_DAL2_Alignment_Viewer.py'
    frame=table();summary,members=summarize_table(frame,REFERENCE,ANNOTATIONS)
    app=AppTest.from_file(str(page))
    app.session_state['dal2_summary']=(summary,members,len(frame))
    app.session_state['dal2_dataset']=(frame,REFERENCE,ANNOTATIONS,True,.7,80,.8)
    app.run(timeout=30)
    assert not app.exception
    assert '5 analyzed rows' in app.success[0].value
    assert app.button[0].label=='Analyze ALL Excel/CSV rows'
    app.selectbox[0].set_value('uORF2 ATG retained').run()
    assert len(app.multiselect[0].options)==1
    next(b for b in app.button if b.label.startswith('Build interactive')).click().run(timeout=30)
    assert not app.exception
    assert len(app.session_state['dal2_view'][2])==1
    assert len(app.session_state['dal2_summary'][1])==5


def test_context_filter_reports_evidence_and_changes_assignment():
    frame=table().iloc[2:3].copy()
    summary,members=summarize_table(frame,REFERENCE,ANNOTATIONS)
    assert members.group_context_identity.iloc[0]==1
    assert members.group_near_start_identity.iloc[0]==1
    assert bool(members.group_placement_consistent.iloc[0])
    assert members.group_context_covered_bases.iloc[0]==50
    # An intact ATG alone is not sufficient if its downstream alignment is weak.
    from nanopore.summary import retained_group
    from nanopore.alignment import align_consensus
    row=align_consensus(REFERENCE[13:563],REFERENCE,ANNOTATIONS)
    row['start_site_support']['uORF2 ATG']['near_identity']=.5
    assert retained_group(row,ANNOTATIONS)!='uORF2 ATG retained'
