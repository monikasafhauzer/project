from streamlit.testing.v1 import AppTest


def test_app_initial_and_validation():
    app=AppTest.from_file(str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    assert not app.exception
    assert app.title[0].value=='Nanopore 5′ RACE Amplicon Analyzer'
    app.button[0].click().run()
    assert 'at least one sample' in app.error[0].value
    app.radio(key='analysis_mode').set_value('Reference matching').run()
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


def test_local_file_mode_runs_and_renders(tmp_path):
    from pathlib import Path
    from test_local import write_fastq,settings
    from test_analysis import random_dna
    from nanopore.io import reverse_complement
    cfg=settings()
    seq=cfg['anchor']+random_dna(500,12)+reverse_complement(cfg['gene'])
    path=tmp_path/'reads.fastq.gz';write_fastq(path,[seq,seq])
    app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py')).run(timeout=30)
    app.radio(key='input_mode').set_value('Read large files from this computer').run()
    assert not app.exception
    app.text_area(key='paths_1').set_value(str(path))
    for widget in app.text_input:
        if 'Folder to save' in widget.label: widget.set_value(str(tmp_path/'results'))
        elif 'anchor primer' in widget.label: widget.set_value(cfg['anchor'])
        elif 'Gene-specific' in widget.label: widget.set_value(cfg['gene'])
    app.button[0].click().run(timeout=30)
    assert not app.exception
    assert not app.error
    from test_jobs import wait_for_job
    job_id=app.session_state['active_job_reads']
    assert wait_for_job(job_id)['state']=='completed'
    # Reopen with a completely fresh browser session; uploads/session state are gone.
    app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py')).run(timeout=30)
    next(b for b in app.button if b.label=='Load completed results').click().run(timeout=30)
    assert not app.exception
    result,metadata=app.session_state['analysis']
    assert result['stats'].input_reads.tolist()==[2]
    assert metadata['local_files']
    assert len(app.dataframe)>=2
