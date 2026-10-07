"""Full-table grouping by first confidently retained annotated initiation site."""
import math
import pandas as pd
from .alignment import align_consensus, validate_annotations, report_frame

UNCERTAIN = 'Uncertain / altered or incompletely covered start'
NO_START = 'No intact annotated start observed'


def retained_group(row, annotations):
    if row['status'] != 'Aligned':
        return UNCERTAIN
    starts = sorted((a for a in annotations if a.get('kind') == 'start'), key=lambda a: a['start'])
    for start in starts:
        state = row['features'].get(start['label'], 'uncertain')
        if state == 'intact':
            return start['label'] + ' retained'
        if state != 'outside observed 5′ sequence':
            # Do not silently skip a possibly retained/mutated earlier start.
            return UNCERTAIN
    return NO_START


def numeric(value, field, count=False):
    if value is None or pd.isna(value):
        return None
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{field} must contain nonnegative numbers or empty cells.') from exc
    if not math.isfinite(number) or number < 0 or (count and not number.is_integer()):
        raise ValueError(f'{field} must contain nonnegative '+('whole numbers.' if count else 'finite numbers.'))
    if not count and number > 100:
        raise ValueError(f'{field} must use percentages from 0 to 100, as exported by the app.')
    return int(number) if count else number


def summarize_table(frame, reference, annotations, trim=True, min_identity=.7, min_span=80, progress=None):
    validate_annotations(reference, annotations)
    seen = set()
    prepared = []
    for record in frame.to_dict('records'):
        metadata = {k: (None if pd.isna(v) else v) for k,v in record.items() if k != 'consensus_dna'}
        sample = str(metadata.get('sample') or 'Unspecified sample')
        name = str(metadata['amplicon'])
        if (sample,name) in seen:
            raise ValueError(f'Duplicate amplicon {name} in {sample}. Remove duplicate rows to avoid double counting.')
        seen.add((sample,name))
        metadata['sample'] = sample
        metadata['supporting_reads'] = numeric(metadata.get('supporting_reads'), 'supporting_reads', True)
        for field in ('percent_quality_filtered','percent_assigned'):
            metadata[field] = numeric(metadata.get(field), field)
        prepared.append((record['consensus_dna'],metadata))
    members=[]
    for index,(sequence,metadata) in enumerate(prepared):
        result=align_consensus(sequence,reference,annotations,trim,min_identity,min_span)
        result.pop('bases',None)
        result.pop('insertions',None)
        row=dict(id=f'row_{index+1}',name=str(metadata['amplicon']),metadata=metadata,**result)
        row['retained_start_group']=retained_group(row,annotations)
        members.append(row)
        if progress: progress((index+1)/len(prepared))
    details=report_frame(members)
    summary=[]
    for sample, sample_rows in details.groupby('sample',sort=False):
        total=sample_rows.supporting_reads.sum() if sample_rows.supporting_reads.notna().all() else None
        for group,group_rows in sample_rows.groupby('retained_start_group',sort=False):
            count=group_rows.supporting_reads.sum() if group_rows.supporting_reads.notna().all() else None
            summary.append(dict(sample=sample,retained_start_group=group,amplicon_clusters=len(group_rows),
                supporting_reads=count,percent_of_imported_supporting_reads=100*count/total if count is not None and total else None,
                percent_quality_filtered_from_input=group_rows.percent_quality_filtered.sum() if group_rows.percent_quality_filtered.notna().all() else None,
                percent_assigned_from_input=group_rows.percent_assigned.sum() if group_rows.percent_assigned.notna().all() else None,
                clusters_with_unaligned_5prime_extension=int(group_rows.get('upstream_sequence',pd.Series(dtype=str)).fillna('').astype(str).ne('').sum()),
                imported_clusters_in_sample=len(sample_rows)))
    return pd.DataFrame(summary),details
