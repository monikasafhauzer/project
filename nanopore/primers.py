from dataclasses import dataclass
import regex
from .io import reverse_complement

IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}


@dataclass
class Amplicon:
    read_id: str
    sequence: str
    status: str
    reversed: bool
    edits: int


def sites(sequence, primer, edits):
    pattern = regex.compile('(?:' + ''.join(IUPAC[x] for x in primer) + '){e<=' + str(edits) + '}', regex.IGNORECASE)
    # Enumerate overlapping positions; fuzzy counts include substitutions and indels.
    return [(m.start(), m.end(), sum(m.fuzzy_counts)) for m in pattern.finditer(sequence, overlapped=True) if m.end()>m.start()]


def extract(read, anchor, gene_primer, edits=2, keep_single=False, full=False):
    # Gene primer is entered in oligo 5′→3′ orientation; its RC occurs downstream.
    downstream = reverse_complement(gene_primer)
    candidates = []
    singles = []
    for reversed_, seq in [(False, read.sequence), (True, reverse_complement(read.sequence))]:
        left, right = sites(seq, anchor, edits), sites(seq, downstream, edits)
        for a in left:
            for b in right:
                if a[1] < b[0]:
                    fragment = seq[a[0]:b[1]] if full else seq[a[1]:b[0]]
                    candidates.append((a[2]+b[2], -(b[0]-a[1]), reversed_, fragment))
        for a in left:
            singles.append((a[2], reversed_, 'anchor_only', seq if full else seq[a[1]:]))
        for b in right:
            singles.append((b[2], reversed_, 'gene_only', seq if full else seq[:b[0]]))
    if candidates:
        error, _, rev, fragment = min(candidates)
        return Amplicon(read.id, fragment, 'both', rev, error)
    if keep_single and singles:
        error, rev, status, fragment = min(singles)
        if fragment:
            return Amplicon(read.id, fragment, status, rev, error)
    return None
