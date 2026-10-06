import gzip
import io
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np


@dataclass
class Read:
    id: str
    sequence: str
    quality: float | None
    source: str


def dna(text: str) -> str:
    sequence = ''.join(text.split()).upper()
    if not sequence or set(sequence) - set('ACGTRYSWKMBDHVN'):
        raise ValueError('DNA must contain only IUPAC DNA letters (A, C, G, T and ambiguity codes).')
    return sequence


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def load_reads(files, sample):
    reads = []
    for file in files:
        name = file.name.lower()
        if not name.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz')):
            raise ValueError(f'{file.name}: unsupported format. Basecall POD5 with Dorado first.')
        raw = file.getvalue()
        if name.endswith('.gz'):
            raw = gzip.decompress(raw)
            name = name[:-3]
        fmt = 'fastq' if name.endswith(('.fastq', '.fq')) else 'fasta'
        try:
            records = list(SeqIO.parse(io.StringIO(raw.decode('utf-8-sig')), fmt))
            if not records:
                raise ValueError('No sequence records found.')
            for record in records:
                sequence = dna(str(record.seq))
                scores = record.letter_annotations.get('phred_quality')
                quality = None if scores is None else float(-10 * np.log10(np.mean(10.0 ** (-np.asarray(scores) / 10))))
                reads.append(Read(f'{sample}_r{len(reads)+1}', sequence, quality, file.name))
        except (ValueError, UnicodeError) as exc:
            raise ValueError(f'{file.name}: invalid {fmt.upper()} input: {exc}') from exc
    return reads


def references(text):
    if not text.strip():
        return {}
    if text.lstrip().startswith('>'):
        result = {}
        for record in SeqIO.parse(io.StringIO(text.strip()), 'fasta'):
            if record.id in result:
                raise ValueError(f'Duplicate reference name: {record.id}')
            result[record.id] = dna(str(record.seq))
        if not result:
            raise ValueError('No reference sequences found.')
        return result
    return {'reference_1': dna(text)}


def fasta(sequences):
    return ''.join(f'>{name}\n{sequence}\n' for name, sequence in sequences)
