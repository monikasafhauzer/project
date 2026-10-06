import gzip
import io
import random
import pytest
from nanopore.io import Read,load_reads,references,reverse_complement
from nanopore.primers import extract
from nanopore.discovery import discover
from nanopore.matching import assign_hits,match_references
from nanopore.pipeline import filter_reads


class Upload(io.BytesIO):
    def __init__(self,data,name):
        super().__init__(data); self.name=name


def random_dna(n,seed):
    r=random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def test_formats_quality_ids():
    files=[Upload(gzip.compress(b'@same\nACGT\n+\nIIII\n'),'a.fastq.gz'),Upload(b'>same\nACGT\n','b.fasta')]
    reads=load_reads(files,'S1')
    assert len(reads)==2 and reads[0].id!=reads[1].id
    assert reads[0].quality==pytest.approx(40) and reads[1].quality is None
    assert len(filter_reads(reads,7,1,10,False))==1
    assert len(filter_reads(reads,7,1,10,True))==2
    with pytest.raises(ValueError): load_reads([Upload(b'@r\nAAAA\n+\nII\n','bad.fastq')],'S')
    with pytest.raises(ValueError): references('>r\nACGT\n>r\nACGT')


def test_primer_orientation_edits_single():
    anchor='ACGTTGCACTGATCGA'; gene='TTCAGGCTAACGTGAC'; insert=random_dna(300,3)
    seq=anchor+insert+reverse_complement(gene)
    for strand in (seq,reverse_complement(seq)):
        a=extract(Read('r',strand,30,'x'),anchor,gene,0)
        assert a.sequence==insert and a.status=='both'
    for mutant in [anchor[:5]+'A'+anchor[6:], anchor[:5]+anchor[6:],anchor[:5]+'T'+anchor[5:]]:
        assert extract(Read('r',mutant+insert+reverse_complement(gene),30,'x'),anchor,gene,1).status=='both'
    assert extract(Read('r',anchor+insert,30,'x'),anchor,gene,0) is None
    assert extract(Read('r',anchor+insert,30,'x'),anchor,gene,0,True).status=='anchor_only'
    assert extract(Read('r',seq,30,'x'),anchor,gene,0,False,True).sequence==seq


def test_vsearch_different_sequences_same_length():
    from nanopore.primers import Amplicon
    a=random_dna(500,10); b=random_dna(500,20)
    reads=[Amplicon(f'r{i}',seq,'both',False,0) for i,seq in enumerate([a,a,a,b,b])]
    table=discover(reads,'S1',10,.95)
    assert sorted(table.supporting_reads)==[2,3]
    assert table.length.tolist()==[500,500]
    assert table.percent_quality_filtered.sum()==50
    assert table.percent_assigned.sum()==100
    assert set(table.consensus_dna)=={a,b}


def hit(read,ref,cov=.99,identity=.95):
    return dict(read_id=read,reference=ref,reference_coverage=cov,identity=identity,strand='+',score=900)


def test_multi_reference_policy_and_partial():
    assigned,ambiguous=assign_hits([hit('a','x'),hit('a','y'),hit('b','x',.6),hit('c','x',.1),hit('d','x',identity=.5)],.85,.5,.95)
    assert [x['read_id'] for x in assigned]==['b']
    assert assigned[0]['match_type']=='partial'
    assert len(ambiguous)==2 and {x['read_id'] for x in ambiguous}=={'a'}


def test_minimap_real_orientation_partial_ambiguous():
    a=random_dna(1400,101); b=random_dna(1400,202)
    refs={'a':a,'duplicate_a':a,'b':b}
    reads=[Read('r1',a,30,'x'),Read('r2',reverse_complement(b),30,'x'),Read('r3',b[200:1100],30,'x'),Read('r4',random_dna(1400,303),30,'x')]
    table,details,counts=match_references(reads,refs,'S1',.9,.5,.95)
    assert counts==dict(unique=2,ambiguous=1,unmatched=1)
    row=table.set_index('reference').loc['b']
    assert row.matching_reads==2 and row.full_length_reads==1 and row.partial_reads==1
    assert row.percent_quality_filtered==50 and row.percent_reference_assigned==100
    assert details[details.read_id=='r2'].iloc[0].strand=='-'
    assert table.matching_reads.sum()==2


def test_three_sample_discovery_counts_are_independent():
    from nanopore.pipeline import analyze
    anchor='ACGTTGCACTGATCGA'; gene='TTCAGGCTAACGTGAC'; seq=anchor+random_dna(500,90)+reverse_complement(gene)
    samples={f'Sample_{i}':[Read(f's{i}r{j}',seq,30,'x') for j in range(i)] for i in range(1,4)}
    settings=dict(mode='Amplicon discovery',min_quality=7,min_length=100,max_length=20000,retain_fasta=True,
        anchor=anchor,gene=gene,edits=0,keep_single=False,full=False,cluster_identity=.9)
    result=analyze(samples,settings,{})
    assert result['amplicons'].supporting_reads.tolist()==[1,2,3]
    assert result['stats'].quality_filtered_reads.tolist()==[1,2,3]
    assert result['amplicons'].percent_quality_filtered.tolist()==[100,100,100]
    assert len(set(result['amplicons'].consensus_dna))==1
