# Payload Shape Characterization


## FANC.Cas9 (200bp amp, 250 reads)

- payloads (read x ref): **231**

| field | n_nonempty | len p50 | len p90 | len max | total elems | uniq ratio p50 | min | max | cur bytes | narrow bytes | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_positions | 231 | 471 | 471 | 471 | 108576 | 0.476 | -93 | 222 | 868,608 | 217,152 | 4.00x |
| all_insertion_positions | 231 | 2 | 2 | 2 | 462 | 1.000 | 92 | 93 | 3,696 | 462 | 8.00x |
| all_insertion_left_positions | 231 | 1 | 1 | 1 | 231 | 1.000 | 92 | 92 | 1,848 | 231 | 8.00x |
| insertion_positions | 231 | 2 | 2 | 2 | 462 | 1.000 | 92 | 93 | 3,696 | 462 | 8.00x |
| insertion_sizes | 231 | 1 | 1 | 1 | 231 | 1.000 | 92 | 248 | 1,848 | 462 | 4.00x |
| all_deletion_positions | 231 | 221 | 221 | 221 | 51051 | 1.000 | 0 | 222 | 408,408 | 102,102 | 4.00x |
| deletion_positions | 231 | 92 | 92 | 92 | 21252 | 1.000 | 0 | 91 | 170,016 | 21,252 | 8.00x |
| deletion_sizes | 231 | 1 | 1 | 1 | 231 | 1.000 | 92 | 92 | 1,848 | 231 | 8.00x |
| all_substitution_positions | 172 | 1 | 1 | 1 | 172 | 1.000 | 92 | 93 | 1,376 | 172 | 8.00x |
| substitution_positions | 172 | 1 | 1 | 1 | 172 | 1.000 | 92 | 93 | 1,376 | 172 | 8.00x |
| substitution_values | 172 | 1 | 1 | 1 | 172 | 1.000 | - | - | 172 | 172 | 1.00x |
| all_substitution_values | 172 | 1 | 1 | 1 | 172 | 1.000 | - | - | 172 | 172 | 1.00x |
| **TOTAL** | | | | | | | | | **1,463,064** | **343,042** | **4.26x** |

**Per-read list-field bytes** (int arrays only + str): median current = 6,346 B, median narrow = 1,488 B (median ratio 4.26x; p90 current = 6,346 B, p90 narrow = 1,488 B)

**ref_positions structure**: insertion-sentinel (negative) fraction p50=0.527 max=0.527; repeat (deletion) fraction p50=0.524 max=0.524. (Low neg+repeat => mostly a clean +1 ramp = highly compressible / reconstructable from aln_ref.)

**all_deletion_positions range expansion**: ratio = len / (2*n_ranges) p50=55.25 max=55.25 (across 231 reads with deletions). Ratio > 1 means range expansion is wasting bytes vs storing (start,end) pairs.


## HEK3.Cas9 (synthetic amp, 250 reads)

- payloads (read x ref): **249**

| field | n_nonempty | len p50 | len p90 | len max | total elems | uniq ratio p50 | min | max | cur bytes | narrow bytes | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_positions | 249 | 468 | 468 | 482 | 115398 | 0.502 | -118 | 233 | 923,184 | 230,796 | 4.00x |
| all_insertion_positions | 249 | 2 | 2 | 2 | 498 | 1.000 | 117 | 118 | 3,984 | 498 | 8.00x |
| all_insertion_left_positions | 249 | 1 | 1 | 1 | 249 | 1.000 | 117 | 117 | 1,992 | 249 | 8.00x |
| insertion_positions | 249 | 2 | 2 | 2 | 498 | 1.000 | 117 | 118 | 3,984 | 498 | 8.00x |
| insertion_sizes | 249 | 1 | 1 | 1 | 249 | 1.000 | 118 | 248 | 1,992 | 498 | 4.00x |
| all_deletion_positions | 249 | 232 | 232 | 232 | 57768 | 1.000 | 0 | 233 | 462,144 | 115,536 | 4.00x |
| deletion_positions | 249 | 117 | 117 | 117 | 29133 | 1.000 | 0 | 116 | 233,064 | 29,133 | 8.00x |
| deletion_sizes | 249 | 1 | 1 | 1 | 249 | 1.000 | 117 | 117 | 1,992 | 249 | 8.00x |
| all_substitution_positions | 159 | 1 | 1 | 1 | 159 | 1.000 | 117 | 118 | 1,272 | 159 | 8.00x |
| substitution_positions | 159 | 1 | 1 | 1 | 159 | 1.000 | 117 | 118 | 1,272 | 159 | 8.00x |
| substitution_values | 159 | 1 | 1 | 1 | 159 | 1.000 | - | - | 159 | 159 | 1.00x |
| all_substitution_values | 159 | 1 | 1 | 1 | 159 | 1.000 | - | - | 159 | 159 | 1.00x |
| **TOTAL** | | | | | | | | | **1,635,198** | **378,093** | **4.32x** |

**Per-read list-field bytes** (int arrays only + str): median current = 6,594 B, median narrow = 1,525 B (median ratio 4.32x; p90 current = 6,610 B, p90 narrow = 1,529 B)

**ref_positions structure**: insertion-sentinel (negative) fraction p50=0.500 max=0.515; repeat (deletion) fraction p50=0.498 max=0.512. (Low neg+repeat => mostly a clean +1 ramp = highly compressible / reconstructable from aln_ref.)

**all_deletion_positions range expansion**: ratio = len / (2*n_ranges) p50=58.00 max=58.00 (across 249 reads with deletions). Ratio > 1 means range expansion is wasting bytes vs storing (start,end) pairs.


## synth_long (2000bp amp, 200 reads)

- payloads (read x ref): **200**

| field | n_nonempty | len p50 | len p90 | len max | total elems | uniq ratio p50 | min | max | cur bytes | narrow bytes | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_positions | 200 | 2004 | 2009 | 2012 | 400909 | 1.000 | -1997 | 1999 | 3,207,272 | 801,818 | 4.00x |
| all_insertion_positions | 197 | 6 | 12 | 16 | 1242 | 1.000 | 1 | 1997 | 9,936 | 2,484 | 4.00x |
| all_insertion_left_positions | 197 | 3 | 6 | 8 | 621 | 1.000 | 1 | 1996 | 4,968 | 1,242 | 4.00x |
| all_deletion_positions | 196 | 4 | 9 | 13 | 912 | 1.000 | 0 | 1999 | 7,296 | 1,824 | 4.00x |
| deletion_positions | 1 | 2 | 2 | 2 | 2 | 1.000 | 1001 | 1002 | 16 | 4 | 4.00x |
| deletion_sizes | 1 | 1 | 1 | 1 | 1 | 1.000 | 2 | 2 | 8 | 1 | 8.00x |
| all_substitution_positions | 70 | 1 | 2 | 3 | 86 | 1.000 | 31 | 1991 | 688 | 172 | 4.00x |
| all_substitution_values | 70 | 1 | 2 | 3 | 86 | 1.000 | - | - | 86 | 86 | 1.00x |
| **TOTAL** | | | | | | | | | **3,230,270** | **807,631** | **4.00x** |

**Per-read list-field bytes** (int arrays only + str): median current = 16,144 B, median narrow = 4,036 B (median ratio 4.00x; p90 current = 16,248 B, p90 narrow = 4,062 B)

**ref_positions structure**: insertion-sentinel (negative) fraction p50=0.002 max=0.006; repeat (deletion) fraction p50=0.000 max=0.003. (Low neg+repeat => mostly a clean +1 ramp = highly compressible / reconstructable from aln_ref.)

**all_deletion_positions range expansion**: ratio = len / (2*n_ranges) p50=0.60 max=2.50 (across 196 reads with deletions). Ratio > 1 means range expansion is wasting bytes vs storing (start,end) pairs.


## synth_long (10000bp amp, 200 reads)

- payloads (read x ref): **200**

| field | n_nonempty | len p50 | len p90 | len max | total elems | uniq ratio p50 | min | max | cur bytes | narrow bytes | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_positions | 200 | 10004 | 10008 | 10012 | 2000873 | 1.000 | -9941 | 9999 | 16,006,984 | 4,001,746 | 4.00x |
| all_insertion_positions | 194 | 6 | 10 | 16 | 1192 | 1.000 | 27 | 9941 | 9,536 | 2,384 | 4.00x |
| all_insertion_left_positions | 194 | 3 | 5 | 8 | 596 | 1.000 | 27 | 9940 | 4,768 | 1,192 | 4.00x |
| all_deletion_positions | 190 | 4 | 9 | 13 | 845 | 1.000 | 26 | 9939 | 6,760 | 1,690 | 4.00x |
| all_substitution_positions | 92 | 1 | 2 | 3 | 116 | 1.000 | 167 | 9988 | 928 | 232 | 4.00x |
| all_substitution_values | 92 | 1 | 2 | 3 | 116 | 1.000 | - | - | 116 | 116 | 1.00x |
| **TOTAL** | | | | | | | | | **16,029,092** | **4,007,360** | **4.00x** |

**Per-read list-field bytes** (int arrays only + str): median current = 80,144 B, median narrow = 20,036 B (median ratio 4.00x; p90 current = 80,224 B, p90 narrow = 20,056 B)

**ref_positions structure**: insertion-sentinel (negative) fraction p50=0.000 max=0.001; repeat (deletion) fraction p50=0.000 max=0.001. (Low neg+repeat => mostly a clean +1 ramp = highly compressible / reconstructable from aln_ref.)

**all_deletion_positions range expansion**: ratio = len / (2*n_ranges) p50=0.50 max=2.50 (across 190 reads with deletions). Ratio > 1 means range expansion is wasting bytes vs storing (start,end) pairs.

