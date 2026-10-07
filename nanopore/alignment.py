"""Reference-anchored consensus comparison, independent of read analysis."""
import io
import pandas as pd
import regex
from Bio.Align import PairwiseAligner
from .io import dna, reverse_complement
from .dal2 import TSO, TSO_TAIL, FORWARD_PRIMER


def load_consensus(name, content):
    if name.lower().endswith('.xlsx'):
        frame=pd.read_excel(io.BytesIO(content),engine='openpyxl')
    elif name.lower().endswith('.csv'):
        frame=pd.read_csv(io.BytesIO(content))
    else:
        raise ValueError('Use a CSV or .xlsx workbook (first worksheet).')
    frame.columns=[str(x).strip() for x in frame.columns]
    missing={'amplicon','consensus_dna'}-set(frame.columns)
    if missing:
        raise ValueError('Missing columns: '+', '.join(sorted(missing)))
    if frame.empty:
        raise ValueError('The table contains no amplicons.')
    if len(frame)>5000:
        raise ValueError('This viewer accepts up to 5,000 rows. Export a selected subset for comparison.')
    for i,value in enumerate(frame.consensus_dna):
        try:
            seq=dna(str(value))
            if pd.isna(value) or len(seq)>20000:
                raise ValueError('Missing sequence or sequence longer than 20,000 bases.')
            frame.loc[frame.index[i],'consensus_dna']=seq
        except ValueError as exc:
            raise ValueError(f'Row {i+2}: {exc}') from exc
    return frame


def validate_annotations(reference, annotations):
    for item in annotations:
        if not item['label'] or not 1<=int(item['start'])<=int(item['end'])<=len(reference):
            raise ValueError('Every annotation needs a label and valid 1-based inclusive reference coordinates.')
        if item.get('kind')=='start' and reference[int(item['start'])-1:int(item['end'])]!='ATG':
            raise ValueError(f"{item['label']}: the annotated start interval is not ATG in this reference.")


def trim_tso(sequence, enabled=True, edits=1):
    if not enabled:
        return sequence, '', 'Adapter trimming disabled; endpoint interpretation may include TSO'
    candidates=[]
    # Only near the physical 5′ boundary; never trim an internal motif.
    for motif in (TSO, FORWARD_PRIMER+TSO_TAIL, TSO_TAIL):
        pattern='(?:'+motif+'){e<='+str(edits)+'}'
        match=regex.search(pattern,sequence[:len(TSO)+8],regex.BESTMATCH)
        if match and match.start()<=5:
            candidates.append((sum(match.fuzzy_counts),-len(motif),match.end()))
    if not candidates:
        return sequence,'','No recognizable TSO boundary; not proof of adapter absence'
    _,_,end=min(candidates)
    return sequence[end:],sequence[:end],'Candidate TSO removed; terminal G ownership/boundary remains uncertain'


def make_aligner():
    aligner=PairwiseAligner(mode='local')
    aligner.match_score=2
    aligner.mismatch_score=-3
    aligner.open_gap_score=-5
    aligner.extend_gap_score=-1
    return aligner


def align_consensus(sequence,reference,annotations,trim=True,min_identity=.7,min_span=80):
    aligner=make_aligner()
    options=[]
    for reverse,oriented in ((False,sequence),(True,reverse_complement(sequence))):
        cleaned,adapter,note=trim_tso(oriented,trim)
        if not cleaned:
            continue
        score=aligner.score(reference,cleaned)
        options.append((score,reverse,cleaned,adapter,note))
    if not options or max(x[0] for x in options)<=0:
        return dict(status='Insufficient alignment',sequence=sequence,features={},bases=[],insertions=[])
    score,reverse,cleaned,adapter,note=max(options,key=lambda x:x[0])
    alignments=iter(aligner.align(reference,cleaned))
    alignment=next(alignments)
    alternate=next(alignments,None)
    ambiguous=alternate is not None and tuple(alternate.coordinates[:,0])!=tuple(alignment.coordinates[:,0])
    coords=alignment.coordinates
    rstart,qstart=map(int,coords[:,0]); rend,qend=map(int,coords[:,-1])
    bases=[]; insertions=[]; matches=columns=0; mapping={}
    for i in range(coords.shape[1]-1):
        ra,qa=map(int,coords[:,i]); rb,qb=map(int,coords[:,i+1])
        if rb>ra and qb>qa:
            for r,q in zip(range(ra,rb),range(qa,qb)):
                base=cleaned[q]; identical=base==reference[r] and base in 'ACGT'
                matches+=identical; columns+=1; mapping[r]=base
                bases.append(dict(position=r+1,base=base,kind='match' if identical else 'mismatch'))
        elif rb>ra:
            columns+=rb-ra
            for r in range(ra,rb):
                mapping[r]='-'; bases.append(dict(position=r+1,base='-',kind='deletion'))
        else:
            columns+=qb-qa
            insertions.append(dict(after=ra,sequence=cleaned[qa:qb]))
    identity=matches/columns if columns else 0
    status='Aligned' if identity>=min_identity and rend-rstart>=min_span else 'Insufficient alignment'
    if status=='Aligned' and ambiguous:
        status='Ambiguous alignment'
    upstream=cleaned[:qstart]; downstream=cleaned[qend:]
    if status!='Aligned':
        endpoint='Uncertain: insufficient reference alignment'
        start=None
    elif upstream and rstart==0:
        endpoint='Extends upstream of reference/uORF1; exact upstream position unknown'
        start=None
    elif upstream:
        endpoint='Uncertain: unaligned 5′ segment before a downstream alignment'
        start=None
    else:
        start=rstart+1
        endpoint='Observed aligned 5′ boundary (not a proven transcription start)'
    features={}
    for feature in annotations:
        a,b=int(feature['start'])-1,int(feature['end'])
        if status!='Aligned':
            state='uncertain'
        elif all(x in mapping for x in range(a,b)):
            observed=''.join(mapping[x] for x in range(a,b))
            state='intact' if observed==reference[a:b] and not any(a<x['after']<b for x in insertions) else 'altered/gapped'
        elif b<=rstart and not upstream:
            state='outside observed 5′ sequence'
        else:
            state='uncertain/not covered'
        features[feature['label']]=state
    retained=[x['label'] for x in annotations if x.get('kind')=='start' and features[x['label']]=='intact']
    return dict(status=status,sequence=cleaned,reverse_complemented=reverse,adapter_removed=adapter,adapter_note=note,
        identity=identity,reference_coverage=(rend-rstart)/len(reference),aligned_reference_start=rstart+1,
        aligned_reference_end=rend,observed_5prime_coordinate=start,endpoint_interpretation=endpoint,
        upstream_sequence=upstream,downstream_sequence=downstream,first_intact_annotated_start=retained[0] if retained else 'None confidently observed',
        features=features,bases=bases,insertions=insertions,alignment_score=score)


def build_rows(frame,reference,annotations,trim=True,min_identity=.7,min_span=80,progress=None):
    validate_annotations(reference,annotations)
    rows=[]
    for index,record in enumerate(frame.to_dict('records')):
        seq=record['consensus_dna']
        result=align_consensus(seq,reference,annotations,trim,min_identity,min_span)
        metadata={key:(None if pd.isna(value) else value) for key,value in record.items() if key!='consensus_dna'}
        rows.append(dict(id=f'row_{index+1}',name=str(record['amplicon']),metadata=metadata,**result))
        if progress: progress((index+1)/len(frame))
    return rows


def report_frame(rows):
    records=[]
    for row in rows:
        report={**row['metadata'],**{k:v for k,v in row.items() if k not in ('metadata','bases','insertions','features')}}
        report.update({f'feature: {name}':state for name,state in row['features'].items()})
        records.append(report)
    return pd.DataFrame(records)
