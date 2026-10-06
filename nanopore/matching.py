from pathlib import Path
from tempfile import TemporaryDirectory
import pandas as pd
from .io import fasta
from .tools import run, threads

SUMMARY_COLUMNS=['sample','reference','reference_length','matching_reads','full_length_reads','partial_reads','percent_quality_filtered','percent_reference_assigned','ambiguous_candidate_reads','assignment_policy']
DETAIL_COLUMNS=['sample','read_id','reference','identity','reference_coverage','strand','match_type','assignment','candidate_references']


def parse_paf(text):
    hits=[]
    for line in text.splitlines():
        f=line.split('\t')
        if len(f)<12:
            raise RuntimeError('Unexpected minimap2 output: malformed PAF record.')
        length=int(f[10])
        hits.append(dict(read_id=f[0],reference=f[5],identity=int(f[9])/length if length else 0,
            reference_coverage=(int(f[8])-int(f[7]))/int(f[6]),strand=f[4],score=int(f[9]),
            ref_start=int(f[7]),ref_end=int(f[8])))
    return hits


def assign_hits(hits, minimum_identity, minimum_coverage, full_coverage):
    grouped={}
    for hit in hits:
        if hit['identity']>=minimum_identity and hit['reference_coverage']>=minimum_coverage:
            key=(hit['read_id'],hit['reference'])
            if key not in grouped or (hit['score'],hit['identity'])>(grouped[key]['score'],grouped[key]['identity']):
                grouped[key]=hit
    reads={}
    for hit in grouped.values():
        reads.setdefault(hit['read_id'],[]).append(hit)
    assigned=[]; ambiguous=[]
    for read_id, candidates in reads.items():
        names='; '.join(sorted(x['reference'] for x in candidates))
        # Conservative policy: any read passing thresholds for >1 reference is excluded.
        for hit in candidates:
            row={**hit,'candidate_references':names,'assignment':'unique' if len(candidates)==1 else 'ambiguous',
                 'match_type':'full_length' if hit['reference_coverage']>=full_coverage else 'partial'}
            (assigned if len(candidates)==1 else ambiguous).append(row)
    return assigned,ambiguous


def match_references(reads, references, sample, identity, coverage, full_coverage=.95):
    hits=[]
    if reads:
        with TemporaryDirectory(prefix='race-map-') as tmp:
            p=Path(tmp)
            (p/'reads.fa').write_text(fasta((r.id,r.sequence) for r in reads))
            (p/'refs.fa').write_text(fasta(references.items()))
            output=run(['minimap2','-x','map-ont','-c','--secondary=yes','-N',str(max(5,len(references))),
                        '-p','0','-t',threads(),str(p/'refs.fa'),str(p/'reads.fa')])
            hits=parse_paf(output)
    assigned,ambiguous=assign_hits(hits,identity,coverage,full_coverage)
    rows=[]
    for name,seq in references.items():
        chosen=[x for x in assigned if x['reference']==name]
        count=len(chosen)
        full=sum(x['match_type']=='full_length' for x in chosen)
        rows.append(dict(sample=sample,reference=name,reference_length=len(seq),matching_reads=count,
            full_length_reads=full,partial_reads=count-full,
            percent_quality_filtered=100*count/len(reads) if reads else 0,
            percent_reference_assigned=100*count/len(assigned) if assigned else 0,
            ambiguous_candidate_reads=sum(x['reference']==name for x in ambiguous),
            assignment_policy='Multiple qualifying references: excluded from unique counts'))
    details=pd.DataFrame([{**x,'sample':sample} for x in assigned+ambiguous],columns=DETAIL_COLUMNS)
    return pd.DataFrame(rows,columns=SUMMARY_COLUMNS),details,dict(unique=len(assigned),
        ambiguous=len({x['read_id'] for x in ambiguous}),unmatched=len(reads)-len(assigned)-len({x['read_id'] for x in ambiguous}))
