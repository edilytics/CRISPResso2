# VCF Parameters Testing and Verification Design

## Overview

Verify that CRISPResso2's `--vcf_output` produces correct VCF files for NHEJ editing outcomes (deletions and insertions at the cut site).

**Strategy**: Golden file comparison using the existing integration test infrastructure in `CRISPResso2_tests`.

**Test flow**:
1. Use `syn-gen` to generate synthetic FASTQ with documented edits
2. Run CRISPResso2 with `--vcf_output --amplicon_coordinates`
3. Compare output VCF against expected VCF using `make <target> test`

**Why golden file comparison**:
- Matches existing test pattern (`make <target> test` diffs against `expected_results/`)
- Catches regressions in VCF format, variant positions, allele frequencies
- Easy to update expected results when intentional changes are made (`make <target> update`)

**Scope**: NHEJ outcomes only (deletions, insertions). Base editing and prime editing deferred to future work.

## Test Cases

### Test 1: `vcf-basic`
Single amplicon, mixed edits.
- ~100 reads with 30% edit rate
- Mix of deletions (various sizes 1-10bp) and insertions (1-5bp)
- Validates basic VCF structure, REF/ALT correctness, position accuracy

### Test 2: `vcf-deletions-only`
Deletion edge cases.
- 100% deletion edits at cut site
- Include edge cases: 1bp deletion, large deletion (20bp+)
- Validates deletion representation in VCF (anchor base + deleted sequence)

### Test 3: `vcf-insertions-only`
Insertion edge cases.
- 100% insertion edits at cut site
- Include 1bp and multi-base insertions
- Validates insertion representation (anchor base, inserted sequence in ALT)

### Test 4: `vcf-multi-amplicon`
Multiple amplicons with different coordinates.
- Two amplicons with different `--amplicon_coordinates` (e.g., "chr1:1000,chr2:5000")
- Verifies CHROM and POS are correctly assigned per amplicon
- Validates that variants from each amplicon appear with correct genomic positions

### Test 5: `vcf-no-edits`
Unedited reads only.
- 0% edit rate - all reads match reference
- Validates that VCF is empty or contains only header (no spurious variants)

## Implementation Workflow

For each test case:

```bash
# 1. Generate synthetic data
cd CRISPResso2_tests/syn-gen
python syn_gen.py --amplicon <SEQ> --guide <GUIDE> \
    --num-reads 100 --edit-rate 0.3 --seed 42 \
    --output-prefix ../cli_integration_tests/inputs/vcf_basic

# 2. Run CRISPResso2 with VCF output
CRISPResso -r1 cli_integration_tests/inputs/vcf_basic.fastq \
    -a <AMPLICON> -g <GUIDE> \
    --vcf_output --amplicon_coordinates "AMPLICON:1" \
    -n vcf-basic --place_report_in_output_folder

# 3. Add test using test_manager.py
python test_manager.py add CRISPResso_on_vcf-basic/
```

`test_manager.py add` will:
- Copy input files to `cli_integration_tests/inputs/`
- Copy results to `cli_integration_tests/expected_results/`
- Add Makefile target automatically
- Update `.gitignore`

After adding all tests:
- Run `make vcf-basic test` to verify diff passes
- Commit the new inputs, expected results, and Makefile changes

## VCF Verification Criteria

**Header correctness**:
- `##fileformat=VCFv4.x` present
- INFO field definitions (AF for allele frequency)
- Contig definitions match amplicon names

**Variant records**:
- CHROM matches amplicon name (or chromosome from `--amplicon_coordinates`)
- POS is 1-based genomic coordinate (amplicon_start + edit_position)
- REF contains anchor base + deleted sequence (for deletions)
- ALT contains anchor base + inserted sequence (for insertions)
- AF (allele frequency) approximately matches expected frequency

**Edge cases caught by diff**:
- Off-by-one position errors
- Incorrect anchor base selection
- Missing variants
- Spurious variants from sequencing errors (should be filtered)

**Manual spot-check** (one-time during initial setup):
- Compare syn-gen's `_edits.tsv` against CRISPResso2's VCF
- Verify a few variants by hand to confirm the expected VCF is correct before committing

## Implementation Steps

1. Generate synthetic test data for each of the 5 test cases using `syn-gen`
2. Run CRISPResso2 with `--vcf_output` for each test case
3. Manually spot-check one VCF against the edits TSV to verify correctness
4. Use `test_manager.py add` for each test to auto-generate Makefile targets
5. Run `make <target> test` for each to confirm diffs pass
6. Commit inputs, expected results, and Makefile changes to the branch

## Future Work

- **Base editing**: Extend `syn-gen` to generate substitution patterns (C→T, A→G), add `vcf-base-editor` test
- **Prime editing**: Extend `syn-gen` for precise small edits, add `vcf-prime-editor` test
- **Multi-amplicon pooled**: Test VCF output from `CRISPRessoPooled` with multiple amplicons
- **bcftools validation**: Optional CI step to run `bcftools stats` on output VCF
