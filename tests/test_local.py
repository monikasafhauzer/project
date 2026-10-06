import gzip
from pathlib import Path
import pandas as pd
import pytest
from nanopore.io import Read, reverse_complement
from nanopore.local import local_path, iter_reads, weighted_median, analyze_local
from nanopore.pipeline import analyze
from test_analysis import random_dna


def settings(mode='Amplicon discovery'):
    return dict(mode=mode,min_quality=7,min_length=100,max_length=20000,retain_fasta=True,
        anchor='ACGTTGCACTGATCGA',gene='TTCAGGCTAACGTGAC',edits=0,keep_single=False,full=False,
        cluster_identity=.9,reference_identity=.9,reference_coverage=.5,full_coverage=.95)


def write_fastq(path, sequences):
    with gzip.open(path,'wt') as handle:
        for i,seq in enumerate(sequences):
            handle.write(f'@r{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')


def test_paths_and_streaming_gzip(tmp_path,monkeypatch):
    assert str(local_path(r'"C:\Users\Name\reads.fastq.gz"'))=='/mnt/c/Users/Name/reads.fastq.gz'
    p=tmp_path/'reads.fastq.gz'
    write_fastq(p,['ACGT']*10)
    monkeypatch.setattr(gzip,'decompress',lambda *args: (_ for _ in ()).throw(AssertionError('whole-file decompression forbidden')))
    reads=iter_reads([p],'S')
    assert iter(reads) is reads
    assert next(reads).id=='S_r1'
    assert len(list(reads))==9
    assert weighted_median({100:3,300:1})==100
    assert weighted_median({100:2,300:2})==200


def test_local_discovery_all_reads_and_disk_reports(tmp_path):
    cfg=settings(); seq=cfg['anchor']+random_dna(500,12)+reverse_complement(cfg['gene'])
    a=tmp_path/'a.fastq.gz'; b=tmp_path/'b.fastq.gz'
    write_fastq(a,[seq,reverse_complement(seq)]); write_fastq(b,[seq])
    samples={'Sample_1':[a,b],'Sample_2':[b]}
    progress=[]
    result,meta=analyze_local(samples,cfg,{},tmp_path/'results',progress.append)
    assert result['stats'].quality_filtered_reads.tolist()==[3,1]
    assert result['amplicons'].supporting_reads.tolist()==[3,1]
    assert result['amplicons'].both_primers.tolist()==[3,1]
    assert result['lengths'].read_count.sum()==4
    root=Path(meta['report_directory'])
    assert (root/'COMPLETE').exists()
    table=pd.read_csv(root/'Sample_1'/'read_statistics.csv')
    assert len(table)==3 and table.read_id.is_unique
    text=(root/'Sample_1'/'preliminary_amplicons.fasta').read_text()
    assert 'Sample_1_amplicon_1|reads=3|preliminary' in text
    assert not (root/'Sample_1'/'extracted.fa').exists()
    # A rerun preserves existing outputs and uses a different folder.
    _, second=analyze_local({'Sample_1':[b]},cfg,{},tmp_path/'results')
    assert second['report_directory']!=meta['report_directory']
    assert len(pd.read_csv(root/'Sample_1'/'read_statistics.csv'))==3


def test_local_reference_batches_preserve_denominators(tmp_path,monkeypatch):
    import nanopore.local as module
    cfg=settings('Reference matching'); a=random_dna(1400,101); b=random_dna(1400,202)
    refs={'a':a,'dup':a,'b':b}
    seqs=[a,reverse_complement(b),b[200:1100],random_dna(1400,303)]
    p=tmp_path/'reads.fastq.gz'; write_fastq(p,seqs)
    original=module.match_references
    batches=[]
    def bounded_batch(reads,*args):
        batches.append(len(reads))
        return original(reads,*args)
    monkeypatch.setattr(module,'BATCH_READS',2)
    monkeypatch.setattr(module,'match_references',bounded_batch)
    result,meta=analyze_local({'Sample_1':[p]},cfg,refs,tmp_path/'results')
    assert batches==[2,2]
    row=result['matches'].set_index('reference').loc['b']
    assert row.matching_reads==2 and row.percent_quality_filtered==50 and row.percent_reference_assigned==100
    stat=result['stats'].iloc[0]
    assert stat.reference_unique_reads==2 and stat.reference_ambiguous_reads==1 and stat.reference_unmatched_reads==1
    detail=pd.read_csv(Path(meta['report_directory'])/'Sample_1'/'match_details.csv')
    assert set(detail.assignment)=={'unique','ambiguous'}


def test_local_truncated_input_is_incomplete(tmp_path):
    p=tmp_path/'bad.fastq.gz'
    with gzip.open(p,'wt') as handle: handle.write('@r\nACGT\n+\nII\n')
    with pytest.raises(ValueError,match='bad.fastq.gz'):
        analyze_local({'Sample_1':[p]},settings(),{},tmp_path/'results')
    runs=list((tmp_path/'results').iterdir())
    assert (runs[0]/'INCOMPLETE').exists() and not (runs[0]/'COMPLETE').exists()


def test_local_fasta_no_quality_or_amplicons(tmp_path):
    p=tmp_path/'reads.fa'; p.write_text('>r\n'+random_dna(500,4)+'\n')
    result,meta=analyze_local({'Sample_1':[p]},settings(),{},tmp_path/'results')
    assert result['stats'].fasta_reads.iloc[0]==1
    assert result['quality_bins'].empty and 'mean_q' in result['quality_bins'].columns
    assert result['amplicons'].empty
