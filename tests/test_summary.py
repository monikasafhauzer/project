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


def test_uncertain_prior_site_not_silently_skipped():
    frame=pd.DataFrame([
        dict(amplicon='divergent_prefix',consensus_dna='ACCTACCGACTACCTA'+REFERENCE[169:563],supporting_reads=10),
        dict(amplicon='mutated_uorf1',consensus_dna=REFERENCE[:1]+'C'+REFERENCE[2:563],supporting_reads=20),
        dict(amplicon='unmatched',consensus_dna='ACGT'*70,supporting_reads=30)])
    summary,members=summarize_table(frame,REFERENCE,ANNOTATIONS,min_identity=.9)
    assert set(members.retained_start_group)=={UNCERTAIN}
    assert summary.supporting_reads.iloc[0]==60


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
