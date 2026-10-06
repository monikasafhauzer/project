"""Sequential local-file workflow: stream gzip, disk reports, bounded mapping batches.

VSEARCH itself still requires memory for the extracted sequences of ONE sample.
"""
import csv
import gzip
import os
import re
import subprocess
import time
from collections import Counter
from pathlib import Path
from tempfile import mkdtemp

import numpy as np
import pandas as pd
from Bio import SeqIO
from .io import Read, dna
from .primers import extract
from .matching import match_references, SUMMARY_COLUMNS, DETAIL_COLUMNS
from .discovery import COLUMNS
from .tools import threads


BATCH_READS = 1000
BATCH_BASES = 5_000_000

SUPPORTED = ('.fastq', '.fq', '.fasta', '.fa', '.fna')


def local_path(text):
    text = text.strip().strip('"')
    if os.name != 'nt' and re.match(r'^[A-Za-z]:[\\/]', text):
        text = '/mnt/' + text[0].lower() + '/' + text[3:].replace('\\', '/')
    return Path(text).expanduser()


def iter_reads(paths, sample):
    number = 0
    for path in paths:
        path = Path(path)
        name = path.name.lower()
        compressed = name.endswith('.gz')
        base = name[:-3] if compressed else name
        if not base.endswith(SUPPORTED):
            raise ValueError(f'{path.name}: use FASTQ/FASTA, optionally gzip compressed. ZIP/POD5 are not supported.')
        fmt = 'fastq' if base.endswith(('.fastq', '.fq')) else 'fasta'
        count = 0
        try:
            opener = gzip.open if compressed else open
            with opener(path, 'rt', encoding='utf-8-sig') as handle:
                for record in SeqIO.parse(handle, fmt):
                    number += 1
                    count += 1
                    sequence = dna(str(record.seq))
                    scores = record.letter_annotations.get('phred_quality')
                    quality = None if scores is None else float(-10 * np.log10(np.mean(10.0 ** (-np.asarray(scores) / 10))))
                    yield Read(f'{sample}_r{number}', sequence, quality, path.name)
            if not count:
                raise ValueError('No sequence records found.')
        except (ValueError, UnicodeError, OSError, EOFError) as exc:
            raise ValueError(f'{path.name}: cannot read {fmt.upper()}: {exc}') from exc


def weighted_median(counts):
    total = sum(counts.values())
    if not total:
        return None
    positions = [(total - 1)//2, total//2]
    values = []
    cumulative = 0
    for length, count in sorted(counts.items()):
        while positions and positions[0] < cumulative + count:
            values.append(length)
            positions.pop(0)
        cumulative += count
        if not positions:
            return sum(values)/len(values)


def cluster_file(input_path, directory, sample, total, assigned, identity, progress):
    """Stream UC and consensus reports; keep only a bounded UI preview."""
    if not assigned:
        return pd.DataFrame(columns=COLUMNS)
    command = ['vsearch', '--cluster_fast', str(input_path), '--id', str(identity), '--iddef', '1',
               '--strand', 'plus', '--minseqlength', '1', '--threads', threads(),
               '--uc', str(directory/'clusters.uc'), '--consout', str(directory/'consensus-native.fa')]
    started = time.monotonic()
    with (directory/'vsearch.log').open('w') as log:
        with subprocess.Popen(command, stdout=log, stderr=log) as process:
            try:
                while True:
                    try:
                        code = process.wait(timeout=30)
                        break
                    except subprocess.TimeoutExpired:
                        if progress:
                            progress(f'{sample}: VSEARCH clustering, {int((time.monotonic()-started)/60)} minutes elapsed…')
            except BaseException:
                process.kill()
                process.wait()
                raise
    if code:
        with (directory/'vsearch.log').open('rb') as log:
            log.seek(max(0, (directory/'vsearch.log').stat().st_size - 4000))
            tail = log.read().decode(errors='replace')
        raise RuntimeError(f'VSEARCH failed. See {directory}/vsearch.log. {tail}')
    support = Counter()
    both = Counter()
    with (directory/'clusters.uc').open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if fields[0] in ('S', 'H'):
                idx = int(fields[1])
                support[idx] += 1
                both[idx] += fields[8].endswith('|both')
    if sum(support.values()) != assigned:
        raise RuntimeError('VSEARCH did not account for every extracted sequence.')
    preview = []
    consensus_count = 0
    with (directory/'amplicons.csv').open('w', newline='') as handle, (directory/'preliminary_amplicons.fasta').open('w') as exported:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        for idx, record in enumerate(SeqIO.parse(directory/'consensus-native.fa', 'fasta')):
            count = support[idx]
            row = dict(sample=sample, amplicon=f'{sample}_amplicon_{idx+1}', length=len(record.seq),
                       supporting_reads=count, percent_quality_filtered=100*count/total,
                       percent_assigned=100*count/assigned, both_primers=both[idx], single_primer=count-both[idx],
                       consensus_dna=str(record.seq), consensus_status='Preliminary VSEARCH consensus; unpolished')
            writer.writerow(row)
            exported.write(f">{row['amplicon']}|reads={count}|preliminary\n{record.seq}\n")
            consensus_count += 1
            if len(preview) < 500:
                preview.append(row)
    if consensus_count != len(support):
        raise RuntimeError('VSEARCH consensus count differs from cluster count.')
    (directory/'consensus-native.fa').unlink()
    return pd.DataFrame(preview, columns=COLUMNS)


def analyze_local(samples, settings, refs, output_root, progress=None):
    """Process samples sequentially, without retaining complete reads or per-read tables."""
    for paths in samples.values():
        for path in paths:
            if not Path(path).is_file():
                raise ValueError(f'File not found: {path}. Paste its full Windows path, including .gz.')
    root = Path(output_root).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    run = Path(mkdtemp(prefix='run-', dir=root))
    stats, lengths, quality_bins, clusters, matches = [], [], [], [], []
    report_paths = []
    try:
        for sample, paths in samples.items():
            directory = run/sample
            directory.mkdir()
            if progress:
                progress(f'{sample}: reading {len(paths)} file(s) sequentially…')
            n = passed = fasta_count = assigned = both_count = qcount = 0
            qsum = 0.0
            hist, qhist, all_lengths = Counter(), Counter(), Counter()
            batch, batch_bases = [], 0
            match_counts = dict(unique=0, ambiguous=0, unmatched=0)
            ref_counts = {name:dict(matching_reads=0, full_length_reads=0, partial_reads=0, ambiguous_candidate_reads=0) for name in refs}
            with (directory/'read_statistics.csv').open('w', newline='') as raw, \
                 (directory/'primer_details.csv').open('w', newline='') as primers, \
                 (directory/'extracted.fa').open('w') as extracted, \
                 (directory/'match_details.csv').open('w', newline='') as detail:
                read_writer = csv.writer(raw)
                read_writer.writerow(['sample','read_id','source_file','length','mean_q','filter_status'])
                primer_writer = csv.writer(primers)
                primer_writer.writerow(['sample','read_id','primer_status','reverse_complemented','primer_edits'])
                detail_writer = csv.DictWriter(detail, fieldnames=DETAIL_COLUMNS)
                detail_writer.writeheader()

                def flush_batch():
                    nonlocal batch, batch_bases
                    if not batch:
                        return
                    summary, details, counts = match_references(batch, refs, sample, settings['reference_identity'], settings['reference_coverage'], settings['full_coverage'])
                    detail_writer.writerows(details.to_dict('records'))
                    for key in match_counts:
                        match_counts[key] += counts[key]
                    for row in summary.to_dict('records'):
                        for key in ref_counts[row['reference']]:
                            ref_counts[row['reference']][key] += row[key]
                    batch, batch_bases = [], 0

                for read in iter_reads(paths, sample):
                    n += 1
                    length = len(read.sequence)
                    all_lengths[length] += 1
                    if read.quality is None:
                        fasta_count += 1
                    else:
                        qcount += 1
                        qsum += read.quality
                        qhist[round(read.quality, 1)] += 1
                    keep = settings['min_length'] <= length <= settings['max_length'] and (settings['retain_fasta'] if read.quality is None else read.quality >= settings['min_quality'])
                    label = 'Quality-filtered' if keep else 'Excluded'
                    hist[(length,label)] += 1
                    read_writer.writerow([sample,read.id,read.source,length,read.quality,label])
                    if keep:
                        passed += 1
                        if settings['mode'] == 'Amplicon discovery':
                            a = extract(read, settings['anchor'], settings['gene'], settings['edits'], settings['keep_single'], settings['full'])
                            primer_writer.writerow([sample,read.id,a.status if a else 'not_retained',a.reversed if a else '',a.edits if a else ''])
                            if a:
                                assigned += 1
                                both_count += a.status == 'both'
                                extracted.write(f'>{a.read_id}|{a.status}\n{a.sequence}\n')
                        else:
                            batch.append(read)
                            batch_bases += length
                            if len(batch) >= BATCH_READS or batch_bases >= BATCH_BASES:
                                flush_batch()
                    if n % 10000 == 0 and progress:
                        progress(f'{sample}: {n:,} reads examined, {passed:,} passed filters…')
                flush_batch()
            stat = dict(sample=sample,input_reads=n,quality_filtered_reads=passed,excluded_reads=n-passed,
                        fasta_reads=fasta_count,median_read_length=weighted_median(all_lengths),mean_read_q=qsum/qcount if qcount else None)
            lengths.extend(dict(sample=sample,length=length,filter_status=label,read_count=count) for (length,label),count in hist.items())
            quality_bins.extend(dict(sample=sample,mean_q=q,read_count=count) for q,count in qhist.items())
            if settings['mode'] == 'Amplicon discovery':
                if progress:
                    progress(f'{sample}: extracted {assigned:,} reads. Clustering this sample only; VSEARCH memory depends on extracted data size…')
                clusters.append(cluster_file(directory/'extracted.fa', directory, sample, passed, assigned, settings['cluster_identity'], progress))
                stat.update(amplicon_assigned_reads=assigned,both_primer_reads=both_count,single_primer_reads=assigned-both_count)
            else:
                rows = []
                for name, sequence in refs.items():
                    counts = ref_counts[name]
                    rows.append(dict(sample=sample,reference=name,reference_length=len(sequence),**counts,
                        percent_quality_filtered=100*counts['matching_reads']/passed if passed else 0,
                        percent_reference_assigned=100*counts['matching_reads']/match_counts['unique'] if match_counts['unique'] else 0,
                        assignment_policy='Multiple qualifying references: excluded from unique counts'))
                table = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
                table.to_csv(directory/'reference_matches.csv',index=False)
                matches.append(table)
                stat.update(reference_unique_reads=match_counts['unique'],reference_ambiguous_reads=match_counts['ambiguous'],reference_unmatched_reads=match_counts['unmatched'])
            stats.append(stat)
            (directory/'extracted.fa').unlink()
            report_paths.extend(str(x) for x in directory.iterdir() if x.suffix in ('.csv','.fasta'))
        combine = lambda frames: pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
        results = dict(stats=pd.DataFrame(stats),lengths=pd.DataFrame(lengths),quality_bins=pd.DataFrame(quality_bins,columns=['sample','mean_q','read_count']),
                       amplicons=combine(clusters),matches=combine(matches),match_details=pd.DataFrame(),primer_details=pd.DataFrame())
        results['stats'].to_csv(run/'sample_statistics.csv',index=False)
        (run/'COMPLETE').write_text('All samples completed successfully.\n')
        return results, dict(report_directory=str(run),report_files=report_paths,local_files=True)
    except BaseException:
        (run/'INCOMPLETE').write_text('This run did not finish. Partial files must not be interpreted as complete results.\n')
        raise
