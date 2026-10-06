from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
import pandas as pd
from .io import fasta
from .tools import run, threads

COLUMNS = ['sample','amplicon','length','supporting_reads','percent_quality_filtered','percent_assigned','both_primers','single_primer','consensus_dna','consensus_status']


def discover(amplicons, sample, filtered_count, identity):
    if not amplicons:
        return pd.DataFrame(columns=COLUMNS)
    by_id = {a.read_id:a for a in amplicons}
    with TemporaryDirectory(prefix='race-cluster-') as tmp:
        p=Path(tmp)
        (p/'reads.fa').write_text(fasta((a.read_id,a.sequence) for a in amplicons))
        run(['vsearch','--cluster_fast',str(p/'reads.fa'),'--id',str(identity),'--iddef','1',
             '--strand','plus','--minseqlength','1','--threads',threads(), '--uc',str(p/'clusters.uc'),
             '--consout',str(p/'consensus.fa')])
        membership={}
        for line in (p/'clusters.uc').read_text().splitlines():
            fields=line.split('\t')
            if fields[0] in ('S','H'):
                membership.setdefault(int(fields[1]),[]).append(fields[8])
        cons=list(SeqIO.parse(p/'consensus.fa','fasta'))
        if len(cons)!=len(membership) or sum(map(len,membership.values()))!=len(amplicons):
            raise RuntimeError('VSEARCH output did not account for all input amplicons; check its output and sequence lengths.')
        rows=[]
        for idx, record in enumerate(cons):
            ids=membership[idx]
            count=len(ids)
            both=sum(by_id[x].status=='both' for x in ids)
            rows.append(dict(sample=sample,amplicon=f'{sample}_amplicon_{idx+1}',length=len(record.seq),
                supporting_reads=count,percent_quality_filtered=100*count/filtered_count,
                percent_assigned=100*count/len(amplicons),both_primers=both,single_primer=count-both,
                consensus_dna=str(record.seq),consensus_status='Preliminary VSEARCH consensus; unpolished'))
        return pd.DataFrame(rows,columns=COLUMNS)
