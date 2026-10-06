from streamlit.testing.v1 import AppTest


def test_app_initial_and_validation():
    app=AppTest.from_file(str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    assert not app.exception
    assert app.title[0].value=='Nanopore 5′ RACE Amplicon Analyzer'
    app.button[0].click().run()
    assert 'at least one sample' in app.error[0].value
    app.radio[0].set_value('Reference matching').run()
    assert not app.exception
    assert app.text_area[0].label == 'Or paste FASTA / a single DNA sequence'


def test_result_rendering_for_both_modes():
    import pandas as pd
    from pathlib import Path
    from nanopore.pipeline import analyze
    from nanopore.io import Read,reverse_complement
    from test_analysis import random_dna
    anchor='ACGTTGCACTGATCGA'; gene='TTCAGGCTAACGTGAC'; seq=random_dna(1400,101)
    settings=dict(mode='Amplicon discovery',min_quality=7,min_length=100,max_length=20000,retain_fasta=True,
        anchor=anchor,gene=gene,edits=0,keep_single=False,full=False,cluster_identity=.9)
    reads={'Sample_1':[Read('r',anchor+seq+reverse_complement(gene),30,'x')]}
    for mode in ('Amplicon discovery','Reference matching'):
        settings['mode']=mode
        settings.update(reference_identity=.85,reference_coverage=.5,full_coverage=.95)
        refs={'target':seq}
        result=analyze(reads,settings,refs)
        app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py'))
        app.session_state['analysis']=(result,dict(settings=dict(settings),references=refs))
        app.run(timeout=30)
        assert not app.exception
        assert len(app.dataframe)>=2
