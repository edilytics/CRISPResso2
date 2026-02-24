# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CRISPResso2 is a bioinformatics pipeline for analyzing CRISPR/Cas9 genome editing outcomes from deep sequencing data. It aligns reads to reference amplicons, quantifies indels/mutations, and generates reports and visualizations.

**Documentation:** https://docs.crispresso.com

## Common Commands

### Installation (Development)
```bash
pip install -e .
```

### Running Tests

Unit tests (pytest):
```bash
pytest tests/unit_tests/

# With coverage
pytest tests --cov CRISPResso2

# Single test file
pytest tests/unit_tests/test_CRISPRessoCORE.py

# Single test
pytest tests/unit_tests/test_CRISPRessoCORE.py::test_function_name
```

Integration tests (from sibling `../CRISPResso2_tests` directory):
```bash
cd ../CRISPResso2_tests

make install       # Installs CRISPResso2 from ../CRISPResso2
make basic         # Run basic test only
make basic test    # Run basic test and diff against expected results
make all           # Run all integration tests
make all test      # Run all tests and diff against expected results
make clean         # Remove test output directories
```

Individual integration test targets:
- `make basic` - Basic core analysis (CRISPResso_on_FANC.Cas9)
- `make params` - Analysis with advanced parameters
- `make batch` - Batch processing
- `make pooled` - Pooled amplicons
- `make wgs` - WGS analysis
- `make compare` - Sample comparison
- `make aggregate` - Aggregate analysis
- `make prime-editor` - Prime editing analysis
- `make base_editor` - Base editing analysis

Add `test` target to diff results against expected: `make basic test`

### Running CRISPResso Tools

```bash
CRISPResso -r1 reads.fastq -a AMPLICON_SEQUENCE -g GUIDE_SEQUENCE
CRISPRessoBatch -bs batch_file.txt -a AMPLICON_SEQUENCE -g GUIDE_SEQUENCE
CRISPRessoPooled -r1 reads.fastq -f amplicons.txt
CRISPRessoWGS -b aligned.bam -r reference.fa -f regions.txt
CRISPRessoCompare sample1_dir/ sample2_dir/
CRISPRessoAggregate -p 'CRISPResso_on_*'
```

## Architecture

### Entry Points (Console Scripts)

Each tool has a corresponding `*CORE.py` module with a `main()` function:

| Command | Module |
|---------|--------|
| `CRISPResso` | `CRISPRessoCORE.py` |
| `CRISPRessoBatch` | `CRISPRessoBatchCORE.py` |
| `CRISPRessoPooled` | `CRISPRessoPooledCORE.py` |
| `CRISPRessoWGS` | `CRISPRessoWGSCORE.py` |
| `CRISPRessoCompare` | `CRISPRessoCompareCORE.py` |
| `CRISPRessoPooledWGSCompare` | `CRISPRessoPooledWGSCompareCORE.py` |
| `CRISPRessoAggregate` | `CRISPRessoAggregateCORE.py` |

### Core Modules

- **`CRISPRessoCORE.py`** (~8,600 lines) - Main analysis engine: read alignment, indel quantification, result aggregation
- **`CRISPRessoShared.py`** - Exception classes, logging utilities, version info, shared helper functions
- **`writers/vcf.py`** - VCF writing, alternate allele mapping, edit processing
- **`CRISPRessoPlot.py`** (~6,000 lines) - All matplotlib/seaborn visualizations
- **`CRISPRessoMultiProcessing.py`** - Parallel processing orchestration

### Cython Modules (Performance-Critical)

- **`CRISPResso2Align.pyx`** - Custom sequence alignment algorithms
- **`CRISPRessoCOREResources.pyx`** - Data structures including `ResultsSlotsDict` for efficient result storage

Pre-compiled `.so` files exist for macOS (x86_64, arm64) and Linux (x86_64). Rebuild with `pip install -e .` if modifying `.pyx` files.

### Report Generation

- **`CRISPRessoReports/CRISPRessoReport.py`** - Jinja2-based HTML report generation
- **Templates:** `CRISPRessoReports/templates/` - HTML templates for each tool type

### Parameter System

- **`args.json`** - Central parameter definitions for all tools. Contains argument names, types, defaults, help text, and which tools each parameter applies to.

## External Dependencies

Required system tools (for pooled/WGS analysis):
- `bowtie2` - Read alignment
- `samtools` - BAM file processing
- `fastp` - Quality filtering (optional)

## Design Documents

See `design_docs/` for detailed write-ups on specific subsystems and past debugging decisions:

- **`LEFT_NORMALIZATION.md`** - VCF indel left-normalization in `writers/vcf.py`: why it's needed, how the fix works, key data structures

## Data Prep Extraction TODOs

Extract non-trivial computation from `CRISPRessoCORE.py` plot section into `CRISPResso2/plots/data_prep.py`, with corresponding `data_func` in `CRISPRessoPro/plots/data_funcs.py` and registration in `builtin_plots.py`. Each is one commit per repo.

### Done
- [x] plot_3a: `prep_indel_size_distribution` — 99th percentile clipping of xmin/xmax
- [x] plot_3b: `prep_frequency_deletions_insertions` — triple 99th percentile clipping (ins/del/mut)
- [x] plot_4a: `prep_amplicon_modifications` — y_max computation, plot_titles
- [x] plot_4b: `prep_modification_frequency` — plot_title with ref_name logic

### Pattern B — Light/medium transforms
- [ ] plot_1d: `prep_dsODN_piechart` — compute labels/sizes from df_alleles for dsODN detection
- [ ] plot_2a: `prep_nucleotide_quilt` — build modification_percentage_summary_df + nuc_df_for_plot from count vectors
- [ ] plot_2b: `prep_nucleotide_quilt_around_sgRNA` — window slicing of 2a's DataFrames + sgRNA interval adjustment (depends on 2a)
- [ ] plot_4e/4f: `prep_global_modifications_reference` — conditional title/root based on ref_name == HDR

### Pattern B — Base editor plots (require `--base_editor_output`)
- [ ] plot_10d: `prep_log_nuc_freqs` — compute plot_quant_window_idxs from include_idxs
- [ ] plot_10e/10f/10g: `prep_conversion_at_sel_nucs` — compute from_nuc_indices, just_sel_nuc_pcts from nuc percentages

### Pattern C — Heavy cross-ref aggregation
- [ ] plot_4g: `prep_hdr_nucleotide_quilt` — build HDR nuc/mod percentage DataFrames across multiple refs
- [ ] plot_5a/6a/8a: `prep_global_frameshift_data` — aggregate frameshift/splice counts across all refs
- [ ] plot_11a: `prep_pe_nucleotide_quilt` — PE nuc quilt DataFrame assembly (same shape as 4g)
- [ ] plot_11b: `prep_pe_nucleotide_quilt_around_sgRNA` — PE window slicing (same shape as 2b)

### Pattern D — Entangled with file writes / heavy DataFrame processing
- [ ] plot_9: `prep_alleles_around_cut` — DataFrame slicing + prep_alleles_table + sgRNA interval recomputation
- [ ] plot_9a: `prep_amino_acid_table` — amino acid DataFrame processing
- [ ] plot_10h: `prep_base_edit_quilt` — base edit DataFrame + prep_alleles_table
- [ ] plot_10i: `prep_base_edit_upset` — alignment + bp_substitutions computation

## Key Constraints

- **Python 3 only**
- **numpy < 2** required (see test_env.yml)
- Cython build requires numpy headers at compile time
