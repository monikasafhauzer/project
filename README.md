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

## Large files and one-file-at-a-time processing

For multi-GB `.fastq.gz` files, select **Read large files from this computer** instead of browser uploads. Paste each full path in the relevant sample box, one path per line. You can supply only Sample 1 to process one file, or supply all three: samples and files are processed sequentially. Windows paths (`C:\Users\Name\Desktop\reads.fastq.gz`) are converted to `/mnt/c/...` automatically under Ubuntu/WSL. The paths refer to the machine running Streamlit, so this mode is intended for trusted local use.

Leave gzip files compressed. Records are read one at a time with Biopython. Filtering, statistics and primer extraction stream to disk. Reference matching uses batches of at most 1,000 reads or roughly 5 million bases (a single unusually long read can exceed the base target); batches are never loaded together. Counts are aggregated across all batches before percentages are calculated. Unique read IDs persist across multiple files in a sample. No downsampling is applied.

Discovery writes extracted sequences to disk, then runs VSEARCH for **one sample at a time**. VSEARCH still loads that sample's extracted sequences into native memory; this is not a fully disk-backed clustering algorithm. Very large or diverse samples can still exceed available RAM or take hours. Windows installed RAM is not necessarily the WSL memory limit. This workflow was validated on synthetic inputs, not the user's multi-GB datasets.

Each run creates a new folder under the selected report directory (default `~/nanopore-results`). Results include complete per-read statistics/primer or alignment CSVs and complete cluster CSV/FASTA reports. Charts use weighted length/Q histograms (all reads, no sampling); Q bins are rounded to 0.1. To keep the browser responsive, discovery previews and their downloadable tables/FASTA contain at most the first 500 clusters per sample. Full results remain on disk. Across-sample plots in this mode cover only those preview clusters. The **Open report folder in Windows Explorer** button opens the complete report directory when running under WSL. Disk space is needed for extracted FASTA, native clustering intermediates and CSV reports; gzip compression ratios vary. Temporary extracted sequences are removed after successful sample analysis.

Progress reports show examined/retained reads every 10,000 records and elapsed VSEARCH time. Runs with errors are marked `INCOMPLETE`; their partial reports must not be treated as finished results. Successful runs have a `COMPLETE` marker. Reruns preserve earlier outputs in separate folders. Local-file provenance records file paths, sizes and modification times rather than hashing entire multi-GB files; upload-mode provenance includes SHA-256 hashes.

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

This app processes files in the local Python process. VSEARCH/minimap2 use private temporary directories, removed after analysis. Input sequences and results stay in Streamlit session memory; CSV/FASTA downloads are user initiated. Remote users of a shared Streamlit instance send their files to its server: deploy only on a trusted machine. Uploaded files are limited by Streamlit's default 200 MB per file; configure `server.maxUploadSize` deliberately if needed. Browser-upload mode is memory-based and is intended for small datasets. Use local-file mode above for large datasets; native VSEARCH clustering still depends on available memory. Native tools use at most four CPU threads. Upload-mode commands and individual local reference-mapping batches have a 15-minute timeout; local VSEARCH clustering has no fixed timeout. Very permissive primer edits or long reads can be expensive. File SHA-256 hashes, settings, reference DNA and native tool versions are included in downloadable provenance; original sequence identifiers are replaced with unique per-sample read IDs, and source filenames are in the read-length CSV.

## Development and validation

```bash
pip install -r requirements-dev.txt
python -m pytest -q
```

Tests exercise compressed FASTQ/FASTA, malformed records, quality filtering, primer substitutions/indels and both orientations, single-primer behavior, equal-length distinct clusters, real VSEARCH consensus, real minimap2 full/partial/reverse-strand/ambiguous mappings and Streamlit startup/input validation. Native tools must be on PATH.

Modules are separated into `io`, `primers`, `discovery`, `matching`, `pipeline` and `tools`. Improved consensus polishing can replace the consensus stage in `discovery.py`; isoform classification should be a separate downstream module with explicit evidence, rather than inferred from lengths or abundance. `app.py` owns only presentation and input validation.

## Annotated DAL2 consensus alignment viewer

Open **DAL2 Alignment Viewer** in the sidebar or use the link on the main page. This is an independent analysis of grouped consensus CSV/XLSX reports, not another read-level run. Required columns: `amplicon`, `consensus_dna`. Existing sample, read-support, percentages and other metadata are preserved. Excel uses the first worksheet. Select up to 100 rows per comparison from a table of up to 5,000 rows; sequences are limited to 20,000 bases, references to 10,000. The default selection is the first 30 rows, visibly listed in the selection control.

The user-supplied 589-base DAL2 reference starts at uORF1. Default ATG sites are 1, 14, 32, 37, 48 (five upstream starts), 100 (main DAL2), 122 and 170 (candidate internal starts). Reverse-primer binding site: 564–589, matching the reverse complement of `GAGTAAACCTTTCTTTAGTAAGGCCG`. Coordinates are reference-local, 1-based, inclusive; they are not genomic positions. The reference and annotation table can be edited together. Start annotations must match ATG. These annotate initiation sites, not experimentally validated full ORFs or translation.

The PCR forward oligo is `CATTGCAAGCAGTGGTATCAAC`. The TSO continues with `GCAGAGTACATrGrGrG`; original primer-only trimming can leave this residual sequence in the exported consensus. The viewer optionally trims a candidate full TSO, PCR-primer-plus-tail or residual tail near the oriented 5′ boundary, allowing one total edit, preserving the removed sequence and interpretation in reports. RNA G residues are represented as DNA G. No internal motif is trimmed. Missing or truncated adapter is not proof of biological sequence; ownership of G residues at the junction remains uncertain. The input consensus file is never edited.

Biopython local alignment searches both orientations, scoring matches +2, mismatches −3, gap openings −5 and extensions −1. Defaults require 70% aligned identity (exact A/C/G/T matches divided by alignment columns, including gaps) and at least 80 aligned reference bases. These thresholds differ from read clustering identity. Equal-score alternate placements detected at different starts are flagged ambiguous. Weak/ambiguous alignments do not receive feature-presence classifications. Unaligned ends are retained; no full-length exact match is required.

A 5′ extension beyond reference base 1 is reported as upstream of the supplied reference, without invented coordinates. A downstream alignment preceded by an unaligned prefix has an uncertain endpoint. An endpoint aligned directly to the reference is an observed sequence boundary, **not a proven transcription start**. Feature states distinguish intact ATG/site, altered/gapped, outside the observed 5′ sequence and uncertain/not covered. Insertions inside an ATG prevent an intact-codon call. Retaining a start codon does not establish a complete uORF or biological isoform.

Drag rows or use ↑/↓ controls; sort by support, aligned start or length. Zoom for DNA letters and hover for bases, insertions and upstream sequences. The amber upstream lane is visually spaced only: prefixes are not aligned to each other and do not have reference coordinates. At most the last 500 upstream bases are drawn; complete extensions are preserved in reports. The order is remembered in browser storage when available. Export arranged SVG, save/restore order JSON, or download a self-contained HTML viewer for offline interaction. These exports follow the viewer's row order. The Streamlit CSV/JSON exports preserve analysis order and full sequence/alignment data. Reference-projection FASTA omits insertions and unaligned ends and must not be interpreted as full amplicon sequences. Imported abundance percentages retain their original denominators; selecting rows does not recalculate abundance.

### Updating the Windows/WSL installation

Wait for any active run to finish; then stop Streamlit with Ctrl+C before updating. From the folder containing `app.py`:

```bash
curl -fL https://raw.githubusercontent.com/monikasafhauzer/project/nanopore-race-app/scripts/update_dal2_viewer.sh -o /tmp/update_dal2_viewer.sh
bash /tmp/update_dal2_viewer.sh
~/.nanopore-app/bin/python -m streamlit run app.py
```

The updater downloads all required files before replacing them, installs Excel support in the prepared environment, and backs up existing affected files under `update-backups/`. Refresh the browser after restarting. Changing pages during an active Streamlit analysis may interrupt that session; finish the read analysis before navigating to the viewer.

### Full-table retained-start summary

Click **Summarize ALL imported amplicons** in the DAL2 viewer to align every imported row independently of the interactive selection. This can take longer than aligning selected rows. Different upstream extensions are combined when uORF1 is confidently intact. Other groups represent the first confidently retained annotated start (uORF2–5, main DAL2, or either internal ATG). An earlier altered, partially covered or uncertain start blocks assignment to a later group; insufficient or ambiguous alignments remain uncertain. These groups describe retained start sites, not exact transcript starts or complete ORFs.

The summary shows amplicon-cluster counts, summed supporting reads, original-input abundance percentages, extension counts, and supporting-read percentages within the imported table, separately per sample. All groups, including uncertainty, contribute to the imported-table denominator. Missing support makes that sample's imported-support denominator unknown; missing percentages remain unknown rather than becoming zero. Duplicate sample/amplicon rows and invalid negative/fractional read counts are rejected to avoid misleading totals. Counts assume the source amplicon clusters have disjoint read memberships.

Expand each group to inspect its members. Download the summary and complete member CSV (including sequences, upstream extensions and feature states). For a complete sample report, import the complete on-disk `amplicons.csv`, not a 500-cluster preview or an Excel subset. The app knows every imported row was analyzed; it cannot establish that the imported file includes every cluster in the original sample. Results remain a labeled snapshot until recomputed. Native read-level analysis and its filtering are unchanged.
