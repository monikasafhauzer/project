import io
import pandas as pd
import pytest
from nanopore.dal2 import REFERENCE,ANNOTATIONS,TSO_TAIL,REVERSE_PRIMER
from nanopore.alignment import align_consensus,load_consensus,build_rows,report_frame,trim_tso
from nanopore.io import reverse_complement
from nanopore.viewer import viewer_html


def test_reference_and_annotated_starts():
    assert len(REFERENCE)==589
    for feature in ANNOTATIONS:
        if feature['kind']=='start': assert REFERENCE[feature['start']-1:feature['end']]=='ATG'
    assert REFERENCE[563:]==reverse_complement(REVERSE_PRIMER)


def test_tso_trim_and_upstream_extension():
    biological='ACCGTTACCGTACCTACGTC'+REFERENCE[:563]
    seq=TSO_TAIL+biological
    r=align_consensus(seq,REFERENCE,ANNOTATIONS)
    assert r['status']=='Aligned'
    assert r['adapter_removed']==TSO_TAIL
    assert r['upstream_sequence']==biological[:20]
    assert r['observed_5prime_coordinate'] is None
    assert 'Extends upstream' in r['endpoint_interpretation']
    assert r['features']['uORF1 ATG']=='intact'
    assert r['features']['Reverse-primer site']=='uncertain/not covered'
    assert r['sequence']==biological


def test_downstream_start_and_both_orientations():
    seq=REFERENCE[169:563]
    for source in (seq,reverse_complement(seq)):
        r=align_consensus(source,REFERENCE,ANNOTATIONS)
        assert r['status']=='Aligned'
        assert r['observed_5prime_coordinate']==170
        assert r['first_intact_annotated_start']=='Internal ATG B'
        assert r['features']['Internal ATG A']=='outside observed 5′ sequence'
        assert r['features']['Internal ATG B']=='intact'


def test_mismatches_and_insertions_preserved():
    seq=list(REFERENCE[:563]);seq[250]='A' if seq[250]!='A' else 'C'
    seq=''.join(seq[:300])+'CCCC'+''.join(seq[300:])
    r=align_consensus(seq,REFERENCE,ANNOTATIONS)
    assert r['status']=='Aligned' and r['identity']<1
    assert any(x['kind']=='mismatch' for x in r['bases'])
    assert ''.join(x['sequence'] for x in r['insertions'])=='CCCC'
    assert r['sequence']==seq
    # An internal TSO-like motif must never cause 5′ trimming.
    internal=REFERENCE[:200]+TSO_TAIL+REFERENCE[200:]
    assert trim_tso(internal)[0]==internal


def test_insufficient_alignment_not_classified():
    r=align_consensus('ACGT'*70,REFERENCE,ANNOTATIONS,min_identity=.99)
    assert r['status']=='Insufficient alignment'
    assert set(r['features'].values())=={'uncertain'}


def test_csv_excel_validation_and_safe_html():
    frame=pd.DataFrame([dict(amplicon='test</script><script>alert(1)</script>',consensus_dna=REFERENCE[:563],supporting_reads=100)])
    loaded=load_consensus('report.csv',frame.to_csv(index=False).encode())
    buf=io.BytesIO();frame.to_excel(buf,index=False)
    assert len(load_consensus('report.xlsx',buf.getvalue()))==1
    rows=build_rows(loaded,REFERENCE,ANNOTATIONS)
    assert report_frame(rows).supporting_reads.iloc[0]==100
    html=viewer_html(REFERENCE,ANNOTATIONS,rows)
    assert 'test</script>' not in html
    assert '\\u003c/script\\u003e' in html
    assert 'ondrop' in html and 'Export arranged SVG' in html
    with pytest.raises(ValueError): load_consensus('bad.csv',b'amplicon,consensus_dna\na,\n')
    bad=[dict(label='wrong',start=2,end=4,kind='start')]
    with pytest.raises(ValueError,match='not ATG'): build_rows(loaded,REFERENCE,bad)


def test_viewer_page_startup():
    from pathlib import Path
    from streamlit.testing.v1 import AppTest
    page=Path(__file__).resolve().parents[1]/'pages'/'1_DAL2_Alignment_Viewer.py'
    app=AppTest.from_file(str(page)).run(timeout=30)
    assert not app.exception
    app.button[0].click().run()
    assert 'Upload' in app.error[0].value
    rows=build_rows(pd.DataFrame([dict(amplicon='test',consensus_dna=REFERENCE[:563])]),REFERENCE,ANNOTATIONS)
    app.session_state['dal2_view']=(REFERENCE,ANNOTATIONS,rows)
    app.run(timeout=30)
    assert not app.exception
    assert len(app.dataframe)==2  # annotation editor plus results table


def test_insertion_inside_atg_is_not_intact():
    seq=REFERENCE[:14]+'CCCC'+REFERENCE[14:563]
    row=align_consensus(seq,REFERENCE,ANNOTATIONS)
    assert row['status']=='Aligned'
    assert row['features']['uORF2 ATG']=='altered/gapped'


def test_repeated_reference_placement_is_ambiguous():
    from test_analysis import random_dna
    seq=random_dna(120,1234)
    row=align_consensus(seq,seq+'TTTTTTTTTT'+seq,[],trim=False)
    assert row['status']=='Ambiguous alignment'
    assert row['observed_5prime_coordinate'] is None
