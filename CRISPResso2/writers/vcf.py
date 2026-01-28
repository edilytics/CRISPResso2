from collections import defaultdict
from contextlib import nullcontext

from CRISPResso2 import CRISPRessoShared


def _ref_length_from_positions(ref_positions):
    """Infer true reference length from the max non-negative position."""
    nonneg = [p for p in ref_positions if p >= 0]
    return (max(nonneg) + 1) if nonneg else 0


def _process_deletions(row, chrom, pos, alt_map):
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    reads = row["#Reads"]

    for (start, end) in row["deletion_coordinates"]:
        left_index = max(1, pos + start - 1)          # absolute 1-based coordinate
        key = (chrom, left_index)
        ref_start = ref_positions.index(start)
        edit_type = "delete"
        try:
            ref_end = ref_positions.index(end)
        except ValueError:  # deletion extends to the end of the sequence, and end is not in ref_positions
            ref_end = ref_positions.index(end - 1)

        if start == 0:  # deletion at the start of a sequence
            edit_type = "delete_start"
            deleted_edit = ref_str[ref_start:ref_end]
            ref_for_key = ref_str[ref_start:ref_end + 1]
        elif ref_end == len(ref_str) - 1:  # deletion at the end of a sequence
            deleted_edit = ref_str[ref_start:ref_end + 1]
            ref_for_key = ref_str[ref_start - 1:ref_end + 1]
        else:  # deletion in the middle of a sequence
            deleted_edit = ref_str[ref_start:ref_end]
            ref_for_key = ref_str[ref_start - 1:ref_end]

        _upsert_edit(
            alt_map,
            key,
            ref_for_key,
            edit_type=edit_type,
            alt_edit=deleted_edit,
            reads=reads,
        )


def _process_insertions(row, chrom, pos, alt_map):
    """Uses ref_positions and insertion_coordinates to find the correct original base."""
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    aln_str = row["Aligned_Sequence"]
    reads = row["#Reads"]

    sizes = row.get("insertion_sizes", []) or []
    coords = row.get("insertion_coordinates", []) or []

    # Infer true reference length in case insertion is at the end of reference
    ref_len = max(p for p in ref_positions if p >= 0) + 1

    for i, (right_anchor_ref_pos, aligned_start) in enumerate(coords):
        left_index = pos + right_anchor_ref_pos
        key = (chrom, left_index)

        # Use the matching size for this insertion
        size = sizes[i] if i < len(sizes) else 0
        ins_edit = aln_str[aligned_start : aligned_start + size]

        # The anchor base for ref_seq: normally at 'right_anchor_ref_pos';
        # for inserts after the last base (right_anchor == ref_len), anchor to the last base.
        if right_anchor_ref_pos == ref_len and ref_len > 0:
            ref_char_idx = ref_positions.index(ref_len - 1)
        else:
            ref_char_idx = ref_positions.index(right_anchor_ref_pos)

        ref_for_key = ref_str[ref_char_idx]  # single-base ref for insertions

        _upsert_edit(
            alt_map,
            key,
            ref_for_key,
            edit_type="insert",
            alt_edit=ins_edit,
            reads=reads,
        )


def _process_substitutions(row, chrom, pos, alt_map):
    """Uses ref_positions and substitution_position to find the correct original base."""
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    aln_str = row["Aligned_Sequence"]
    reads = row["#Reads"]

    for sub_ref_pos in row["substitution_positions"]:
        left_index = pos + sub_ref_pos
        key = (chrom, left_index)

        char_idx = ref_positions.index(sub_ref_pos)
        alt_base = aln_str[char_idx]
        ref_base = ref_str[char_idx]

        _upsert_edit(
            alt_map,
            key,
            ref_base,
            edit_type="sub",
            alt_edit=alt_base,
            reads=reads,
        )


def _upsert_edit(alt_map, key, ref_seq_for_key, edit_type, alt_edit, reads):
    """Helper function to upsert an edit into alt_map[key], always keeping the longest ref_seq seen at this key."""
    entry = alt_map.get(key)
    if entry is None:
        alt_map[key] = {"ref_seq": ref_seq_for_key, "alt_edits": [[edit_type, alt_edit, reads]]}
        return

    # Merge rule
    if edit_type == "delete" or edit_type == 'delete_start':
        target_len = len(alt_edit)
        for existing in entry["alt_edits"]:
            if existing[0] == edit_type and len(existing[1]) == target_len:
                existing[2] += reads
                break
        else:
            entry["alt_edits"].append([edit_type, alt_edit, reads])
    else:
        for existing in entry["alt_edits"]:
            if existing[0] == edit_type and existing[1] == alt_edit:
                existing[2] += reads
                break
        else:
            entry["alt_edits"].append([edit_type, alt_edit, reads])

    # Keep longest ref_seq
    if len(ref_seq_for_key) > len(entry["ref_seq"]):
        entry["ref_seq"] = ref_seq_for_key


def build_alt_map(df_alleles, amplicon_positions):
    """
    Build a nested alt_map keyed by (chrom, left_index).

    For each allele row:
      - Skip unmodified reads (no insertions, deletions, or substitutions).
      - Use amplicon_positions[Reference_Name] -> (chrom, pos), where pos is 1-based
        absolute coordinate of reference index 0.
      - Emit entries:
          key: (chrom, left_index)
          value: {
              "ref_seq": str,      # longest reference context seen at the key
              "alt_edits": [ [edit_type, alt_edit, reads], ... ]
          }

    Merge rules at a given key:
      - deletions: merge entries with the same deleted length (reads sum)
      - insertions: merge entries with the same inserted string
      - substitutions: merge entries with the same alt base

    Edge cases handled:
      - deletion starting at ref index 0 -> no left flank in ref_seq
      - deletion ending at len(ref) -> safe bound computation (no .index(end))
      - multiple insertions per allele -> each uses its own insertion size
      - substitution indexing via ref_positions (robust to hyphen padding)
    """
    alt_map = {}

    for _, row in df_alleles.iterrows():
        if (
            row.get("n_inserted", 0) == 0
            and row.get("n_deleted", 0) == 0
            and row.get("n_mutated", 0) == 0
        ):
            continue

        ref_name = row["Reference_Name"]
        chrom, pos = amplicon_positions[ref_name]

        # Process each edit class, make sure that you always process deletions first, as it will set the proper reference length
        _process_deletions(row, chrom, pos, alt_map)
        _process_insertions(row, chrom, pos, alt_map)
        _process_substitutions(row, chrom, pos, alt_map)

    return alt_map


def _alt_seq_from_edit(ref_seq, edit_type, alt_edit):
    """Build each edit's ALT string using this position's final ref_seq.

    - Insertions: keep left anchor (ref_seq[0]), inject inserted bases, then append any trailing ref context.
        (appending additional context will only occur if there's also a deletion at this position)
    - Deletions: drop exactly len(alt_edit) bases immediately after the left anchor then append any remaining sequences in ref_seq.
        (this will only occur if there's another deletion with a longer span).
    - Substitutions: replace the left anchor base with the alt base, append trailing ref context.

    Parameters
    ----------
    ref_seq : str
        The reference sequence for this key, always including the left anchor base.
        The length of ref_seq is always 1 + the longest deleted span at this position.
        If there are no deletions, ref_seq is a single base (the left anchor).
    edit_type : str
        The type of edit, one of "insert", "delete", "delete_start", or "sub".
    alt_edit : str
        The inserted bases (for "insert"), the deleted bases (for "delete" and "delete_start"),
        or the alt base (for "sub").


    Returns
    -------
    str
        The final alt sequence for this edit.

    """

    alt_len = len(alt_edit)
    ref_len = len(ref_seq)

    if edit_type == "insert":
        # ALT = left_anchor + inserted + trailing_ref
        return ref_seq[0] + (alt_edit or "") + ref_seq[1:]

    if edit_type == "delete":
        if alt_len >= ref_len:
            raise CRISPRessoShared.BadParameterException(
                f"Deletion alt_edit ({alt_edit}) must be at least 1 base pair shorter than ref_seq: ({ref_seq})."
            )
        return ref_seq[0] + ref_seq[alt_len + 1:]

    if edit_type == "delete_start":
        if len(ref_seq) > 1:
            return ref_seq[alt_len:]
        elif len(ref_seq) == 0 or len(ref_seq) == 1:
            raise CRISPRessoShared.BadParameterException(
                f"The ref_seq for deletion at the start of a seq should always have at least 2 basepairs. {ref_seq}."
            )

    if edit_type == "sub":
        if alt_len != 1:
            raise CRISPRessoShared.BadParameterException("Substitution ALT must be a single base.")
        # ALT = alt_base + trailing_ref
        return alt_edit + ref_seq[1:]

    raise CRISPRessoShared.BadParameterException(f"Unknown edit type: {edit_type!r}")


def _write_vcf_line(chrom, pos, ref_seq, alt_edits, num_reads):
    """Takes one (chrom, pos) and returns the relevant VCF record line."""
    if num_reads is None or num_reads <= 0:
        raise CRISPRessoShared.BadParameterException("Error counting total number of reads.")

    alt_counts = defaultdict(int)
    for edit_type, payload, reads in alt_edits:
        alt_seq = _alt_seq_from_edit(ref_seq, edit_type, payload)
        alt_counts[alt_seq] += int(reads)

    # Deterministic ordering: shortest first, then lexicographic
    sorted_alt_seqs = sorted(alt_counts.keys(), key=lambda s: (len(s), s))

    # Allele frequencies
    denom = float(num_reads)
    afs = [f"{alt_counts[alt_seq] / denom:.3f}" for alt_seq in sorted_alt_seqs]

    # Core VCF columns
    info = f"AF={','.join(afs)}"
    record = f"{chrom}\t{pos}\t.\t{ref_seq}\t{','.join(sorted_alt_seqs)}\t.\tPASS\t{info}"

    return record


def vcf_lines_from_alt_map(alt_map, num_reads, amplicon_lens, vcf_path):
    """After being passed a single location, this generates the final alt_seqs and AFs, and returns the single VCF line."""

    def _key_sort(k):
        chrom, pos = k
        try:
            c_key = (0, int(chrom))         # numeric chroms first (e.g., "1", 2)
        except (TypeError, ValueError):
            c_key = (1, str(chrom))         # non‑numeric chroms next, lexicographic (e.g., "chr1", "chr10")
        return (c_key, int(pos))

    counter = 0
    with open(vcf_path, "w") as f:
        f.write('##fileformat=VCFv4.5\n')
        f.write('##source=CRISPResso2\n')
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        for amplicon_name, amplicon_len in amplicon_lens.items():
            f.write(f'##contig=<ID={amplicon_name},length={amplicon_len}>\n')
        base = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        f.write(base + "\n")
        for (chrom, pos) in sorted(alt_map.keys(), key=_key_sort):
            entry = alt_map[(chrom, pos)]
            ref_seq = entry["ref_seq"]
            alt_edits = entry["alt_edits"]
            f.write(_write_vcf_line(chrom, pos, ref_seq, alt_edits, num_reads) + "\n")
            counter += 1
    return counter


def write_vcf_file(df_alleles, ref_names, ref_lens, args, vcf_path):
    """Orchestrates: parse amplicon coordinates, build alt_map, make VCF, and write to file.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        The alleles DataFrame from CRISPResso2 analysis.
    ref_names : list of str
        The list of reference names (amplicon names).
    ref_lens : dict of key str and value int
        A dict where the key is the ref_name and the value is the length
        of the amplicon sequence.
    args : argparse.Namespace
        The command-line arguments, used to get args.amplicon_coordinates.
    vcf_path : str
        The path to write the VCF file to.

    Returns
    -------
    int
        Solely for testing; the number of lines written to the VCF file (including header).
    """
    try:
        all_coords = args.amplicon_coordinates.strip().split(',')
        if len(all_coords) != len(ref_names):
            raise CRISPRessoShared.BadParameterException(
                f"Number of --amplicon_coordinates ({len(all_coords)}) does not match number of amplicons ({len(ref_names)})."
            )

        amplicon_positions = {}
        for i, coord in enumerate(all_coords):
            chrom, pos_str = coord.strip().split(':')
            pos = int(pos_str)
            amplicon_positions[ref_names[i]] = (chrom, pos)

    except ValueError:
        raise CRISPRessoShared.BadParameterException("Invalid format for --amplicon_coordinates.")

    amplicon_lens = {}
    for ref_name, (chrom, pos) in amplicon_positions.items():
        amplicon_lens[chrom] = ref_lens[ref_name]

    # Denominator for AFs
    num_reads = df_alleles['#Reads'].sum()

    alt_map = build_alt_map(df_alleles, amplicon_positions)
    count = vcf_lines_from_alt_map(alt_map, num_reads, amplicon_lens, vcf_path)

    return count
