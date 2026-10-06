import pandas as pd
from .primers import extract
from .discovery import discover
from .matching import match_references


def filter_reads(reads, min_quality, min_length, max_length, retain_fasta):
    return [r for r in reads if min_length<=len(r.sequence)<=max_length and
            (retain_fasta if r.quality is None else r.quality>=min_quality)]


def analyze(samples, settings, refs, progress=None):
    stats=[]; lengths=[]; clusters=[]; matches=[]; details=[]; primer_rows=[]
    for sample, reads in samples.items():
        if progress: progress(f'Analyzing {sample} ({len(reads):,} reads)…')
        filtered=filter_reads(reads,settings['min_quality'],settings['min_length'],settings['max_length'],settings['retain_fasta'])
        retained={r.id for r in filtered}
        for r in reads:
            lengths.append(dict(sample=sample,read_id=r.id,source_file=r.source,length=len(r.sequence),mean_q=r.quality,
                                filter_status='Quality-filtered' if r.id in retained else 'Excluded'))
        q=[r.quality for r in reads if r.quality is not None]
        stat=dict(sample=sample,input_reads=len(reads),quality_filtered_reads=len(filtered),excluded_reads=len(reads)-len(filtered),
                  fasta_reads=sum(r.quality is None for r in reads),median_read_length=pd.Series([len(r.sequence) for r in reads]).median(),
                  mean_read_q=sum(q)/len(q) if q else None)
        if settings['mode']=='Amplicon discovery':
            found=[]
            for read in filtered:
                a=extract(read,settings['anchor'],settings['gene'],settings['edits'],settings['keep_single'],settings['full'])
                primer_rows.append(dict(sample=sample,read_id=read.id,primer_status=a.status if a else 'not_retained',
                    reverse_complemented=a.reversed if a else None,primer_edits=a.edits if a else None))
                if a: found.append(a)
            clusters.append(discover(found,sample,len(filtered),settings['cluster_identity']))
            stat.update(amplicon_assigned_reads=len(found),both_primer_reads=sum(a.status=='both' for a in found),
                        single_primer_reads=sum(a.status!='both' for a in found))
        else:
            summary,detail,counts=match_references(filtered,refs,sample,settings['reference_identity'],settings['reference_coverage'],settings['full_coverage'])
            matches.append(summary); details.append(detail)
            stat.update(reference_unique_reads=counts['unique'],reference_ambiguous_reads=counts['ambiguous'],reference_unmatched_reads=counts['unmatched'])
        stats.append(stat)
    combine=lambda frames: pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    return dict(stats=pd.DataFrame(stats),lengths=pd.DataFrame(lengths),amplicons=combine(clusters),
                matches=combine(matches),match_details=combine(details),primer_details=pd.DataFrame(primer_rows))
