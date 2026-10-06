# Nanopore 5′ RACE Amplicon Analyzer

A local Streamlit application for independently analyzing up to three Oxford Nanopore 5′ RACE PCR samples, including multiple sequencing files and multiple amplicons per sample. Python 3.11 or newer is required. The original biology exercise remains in `legacy/code.py` and `legacy/test_data.csv`; it is unrelated to this application.

## Install and start

Recommended on Linux/macOS: install the native tools with Conda/Miniforge (Bioconda), then the Python application dependencies:

```bash
conda create -n nanopore-race -c conda-forge -c bioconda python=3.11 minimap2 vsearch pip
conda activate nanopore-race
pip install -r requirements.txt
streamlit run app.py
```

Alternatively, install Python 3.11+, minimap2 and VSEARCH with your system package manager, then:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
minimap2 --version
vsearch --version
streamlit run app.py
```

On Windows use WSL2 with Linux and the same instructions. Both external executables must be on the PATH of the Streamlit process. The sidebar reports tool versions and actionable missing-dependency errors. Discovery requires VSEARCH; reference matching requires minimap2. No API credentials or remote analysis services are needed. Open the address printed by Streamlit on your own computer.

### Prepared cloud machine

Python dependencies are installed in `/workspace/.onboarding/nanopore-venv`; native tools are in `/workspace/.onboarding/nanopore-tools/bin`. From `/workspace/project`:

```bash
export PATH=/workspace/.onboarding/nanopore-tools/bin:$PATH
/workspace/.onboarding/nanopore-venv/bin/python -m streamlit run app.py --server.headless true
```

## A molecular biologist's workflow

1. Upload FASTQ, FASTQ.gz or FASTA files into Sample 1, 2 and 3. Multiple files for one sample are pooled only within that sample. At least one sample is required. FASTA.gz, .fq, .fa and .fna are also supported. **POD5 files must first be basecalled with Dorado**; this app does not run basecalling.
2. Choose discovery or reference matching. Set read length and minimum Q score. FASTQ mean Q is calculated from mean base error probability, not the arithmetic mean of Phred scores. FASTA quality is unknown: explicitly include or exclude those reads with the checkbox.
3. For discovery, enter both actual primer oligos in their **5′→3′ direction**. The anchor occurs upstream and the reverse complement of the gene-specific primer occurs downstream on normalized reads. Both strands are searched. The edit budget per primer includes substitutions, insertions and deletions; IUPAC primer codes are supported. At least eight bases per primer are required. Both-primer hits are prioritized, then lowest total edits, then longest interior. Repeated internal primer sites can create ambiguous boundaries: inspect resulting sequences before interpretation.
4. By default only both-primer reads are used and primer sites are removed. Optionally keep single-primer flanks or primer-containing amplicons. Single-primer reads lack a validated second boundary. Counts of both- and single-primer reads are reported separately.
5. For reference matching, upload FASTA or paste DNA/FASTA. Multiple uniquely named references are supported. Uploaded and pasted references are combined. Adjust identity, minimum coverage and full-length coverage threshold.
6. Click **Analyze samples**, review plots and tables, and download CSV, FASTA and provenance/settings JSON. Results remain a labeled snapshot until you rerun; editing controls does not change existing results.

## What the analysis means

### Discovery

Each sample is filtered independently, reads are normalized by primer direction, and extracted sequences are clustered using VSEARCH `--cluster_fast --iddef 1 --strand plus`. Clusters depend on sequence identity across the global alignment, including terminal gaps, **never length alone**. Clustering is greedy and is not biological isoform classification. VSEARCH `--consout` supplies a **preliminary, unpolished consensus** per cluster; it is not a medaka/Racon-polished reconstruction. Single-primer fragments may cluster separately from complete amplicons. No chimera detection is performed.

Reports contain consensus DNA, length, supporting reads, both-/single-primer support and two explicitly different percentages:

- `percent_quality_filtered`: supporting reads / all quality-filtered reads in that sample.
- `percent_assigned`: supporting reads / all reads retained and assigned to amplicon clusters in that sample.

Cross-sample abundance plots group only **exact consensus DNA equality**. Sequencing errors may separate related biological amplicons. Different sequences of equal length are kept separate.

### Reference matching

Minimap2 uses `-x map-ont -c --secondary=yes -p 0`, with enough secondary slots for the supplied references, searching both read orientations. Identity = matching bases / alignment block length, including gaps. Reference coverage = aligned target span / reference length. For multiple alignments to the same reference, the alignment with most matching bases is used; segments are not stitched. A read is a full-length reference match only if it meets the user-set full-coverage threshold (default 95%); other qualifying matches are partial. This means nearly full reference coverage, not proof of a full transcript.

**Multi-reference reads are not silently double counted.** If more than one reference passes both thresholds, the read is ambiguous and excluded from unique assignment totals. Candidate references remain in the per-read CSV. Candidate ambiguity columns can overlap and must not be summed as unique reads. Identical references intentionally produce ambiguous assignments. Minimap2 is heuristic: very short, repetitive or divergent sequences can be missed; these are unmatched rather than proven absent.

- `percent_quality_filtered`: unique reads matching a reference / all quality-filtered reads.
- `percent_reference_assigned`: unique reads matching a reference / all uniquely reference-assigned reads.
- Unique, ambiguous and unmatched counts reconcile to the quality-filtered total per sample.

**PCR-derived read abundance is not original RNA molecule abundance.** Amplification bias, sequencing errors and chimeras can affect both modes. Do not call clusters, preliminary consensus sequences or artifacts genuine biological isoforms without orthogonal validation.

## Resources, privacy and limitations

This app processes files in the local Python process. VSEARCH/minimap2 use private temporary directories, removed after analysis. Input sequences and results stay in Streamlit session memory; CSV/FASTA downloads are user initiated. Remote users of a shared Streamlit instance send their files to its server: deploy only on a trusted machine. Uploaded files are limited by Streamlit's default 200 MB per file; configure `server.maxUploadSize` deliberately if needed. Parsing and tables are memory-based, not a production streaming pipeline. Large datasets should be filtered/subsampled first. Native tools use at most four CPU threads and a 15-minute command timeout. Very permissive primer edits or long reads can be expensive. File SHA-256 hashes, settings, reference DNA and native tool versions are included in downloadable provenance; original sequence identifiers are replaced with unique per-sample read IDs, and source filenames are in the read-length CSV.

## Development and validation

```bash
pip install -r requirements-dev.txt
python -m pytest -q
```

Tests exercise compressed FASTQ/FASTA, malformed records, quality filtering, primer substitutions/indels and both orientations, single-primer behavior, equal-length distinct clusters, real VSEARCH consensus, real minimap2 full/partial/reverse-strand/ambiguous mappings and Streamlit startup/input validation. Native tools must be on PATH.

Modules are separated into `io`, `primers`, `discovery`, `matching`, `pipeline` and `tools`. Improved consensus polishing can replace the consensus stage in `discovery.py`; isoform classification should be a separate downstream module with explicit evidence, rather than inferred from lengths or abundance. `app.py` owns only presentation and input validation.
