"""Full-table grouping by first confidently retained annotated initiation site."""
import math
import pandas as pd
from .alignment import align_consensus, validate_annotations, report_frame

UNCERTAIN = 'Uncertain / altered or incompletely covered start'
NO_START = 'No intact annotated start observed'


def retained_group(row, annotations, context_identity=.8):
    if row['status'] not in ('Aligned','Ambiguous alignment'):
        return UNCERTAIN
    starts = sorted((a for a in annotations if a.get('kind') == 'start'), key=lambda a: a['start'])
    unresolved=False
    for start in starts:
        # A non-reference prefix does not make an earlier reference uORF retained.
        if int(start['end']) < row.get('aligned_reference_start', 1):
            continue
        state = row['features'].get(start['label'], 'uncertain')
        context = row.get('start_site_support', {}).get(start['label'], {})
        context_good = context.get('covered_bases',0) >= min(40,context.get('window_bases',50)) and context.get('identity',0) >= context_identity and context.get('near_covered_bases',0)>=10 and context.get('near_identity',0)>=max(.85,context_identity)
        if state in ('intact','altered/gapped') and context_good and not context.get('placement_consistent',False):
            return UNCERTAIN
        supported = context_good and context.get('placement_consistent',False)
        if state == 'intact' and supported:
            return start['label'] + ' retained'
        if state == 'intact' and not supported:
            unresolved=True
        if state == 'altered/gapped' and supported:
            # A well-aligned but mutated earlier ATG may be a consensus error.
            return UNCERTAIN
        if state == 'uncertain/not covered' and int(start['start']) < row.get('aligned_reference_start',1) <= int(start['end']) and not row.get('upstream_sequence'):
            return UNCERTAIN
    return UNCERTAIN if unresolved else NO_START


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


def summarize_table(frame, reference, annotations, trim=True, min_identity=.7, min_span=80, progress=None, context_identity=.8):
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
        group=retained_group(result,annotations,context_identity)
        result.pop('bases',None)
        result.pop('insertions',None)
        row=dict(id=f'row_{index+1}',name=str(metadata['amplicon']),metadata=metadata,**result)
        row['retained_start_group']=group
        row['grouping_minimum_context_identity']=context_identity
        chosen=group.removesuffix(' retained')
        evidence=result.get('start_site_support',{}).get(chosen,{})
        row['group_context_identity']=evidence.get('identity')
        row['group_near_start_identity']=evidence.get('near_identity')
        row['group_context_covered_bases']=evidence.get('covered_bases')
        row['group_placement_consistent']=evidence.get('placement_consistent')
        row['grouping_basis']='First intact annotated ATG with supported downstream reference context; non-reference prefix does not establish earlier uORFs'
        row['upstream_interpretation']='Different/unaligned 5′ sequence; origin and exact endpoint unresolved' if row.get('upstream_sequence') else 'No unaligned 5′ prefix detected'
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
    group_order={a['label']+' retained':i for i,a in enumerate(sorted((x for x in annotations if x.get('kind')=='start'),key=lambda x:x['start']))}
    group_order[NO_START]=len(group_order);group_order[UNCERTAIN]=len(group_order)
    result=pd.DataFrame(summary)
    if not result.empty:
        result['_group_order']=result.retained_start_group.map(group_order)
        result=result.sort_values(['sample','_group_order'],kind='stable').drop(columns='_group_order').reset_index(drop=True)
    if int(result.amplicon_clusters.sum()) != len(frame):
        raise RuntimeError('Summary did not account for every imported row.')
    return result,details
