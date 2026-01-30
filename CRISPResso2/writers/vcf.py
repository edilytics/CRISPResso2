from collections import defaultdict

from CRISPResso2 import CRISPRessoShared


def _edits_from_deletions(row, chrom, pos):
    """Yield (chrom, vcf_pos, ref, alt, reads) for each deletion in the row.

    VCF deletion representation:
      - Start of sequence (pos 0):  REF = deleted + following_base, ALT = following_base
      - Middle or end of sequence:  REF = anchor + deleted, ALT = anchor
    """
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    reads = row["#Reads"]

    for (start, end) in row["deletion_coordinates"]:
        left_index = max(1, pos + start - 1)
        ref_start = ref_positions.index(start)
        try:
            ref_end = ref_positions.index(end)
        except ValueError:  # deletion extends to the end of the sequence
            ref_end = ref_positions.index(end - 1)

        if start == 0:
            ref_seq = ref_str[ref_start:ref_end + 1]
            alt_seq = ref_seq[-1]
        elif ref_end == len(ref_str) - 1:
            ref_seq = ref_str[ref_start - 1:ref_end + 1]
            alt_seq = ref_seq[0]
        else:
            ref_seq = ref_str[ref_start - 1:ref_end]
            alt_seq = ref_seq[0]

        yield (chrom, left_index, ref_seq, alt_seq, reads)


def _edits_from_insertions(row, chrom, pos):
    """Yield (chrom, vcf_pos, ref, alt, reads) for each insertion in the row.

    VCF insertion: REF = anchor_base, ALT = anchor_base + inserted_bases
    """
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    aln_str = row["Aligned_Sequence"]
    reads = row["#Reads"]

    sizes = row.get("insertion_sizes", []) or []
    coords = row.get("insertion_coordinates", []) or []

    ref_len = max(p for p in ref_positions if p >= 0) + 1

    for i, (right_anchor_ref_pos, aligned_start) in enumerate(coords):
        left_index = pos + right_anchor_ref_pos

        size = sizes[i] if i < len(sizes) else 0
        ins_bases = aln_str[aligned_start:aligned_start + size]

        # Anchor base: normally at right_anchor_ref_pos;
        # for inserts after the last base, anchor to the last base.
        if right_anchor_ref_pos == ref_len and ref_len > 0:
            ref_char_idx = ref_positions.index(ref_len - 1)
        else:
            ref_char_idx = ref_positions.index(right_anchor_ref_pos)

        anchor = ref_str[ref_char_idx]
        yield (chrom, left_index, anchor, anchor + ins_bases, reads)


def _edits_from_substitutions(row, chrom, pos):
    """Yield (chrom, vcf_pos, ref, alt, reads) for each substitution in the row.

    VCF substitution: REF = original_base, ALT = alt_base
    """
    ref_positions = row["ref_positions"]
    ref_str = row["Reference_Sequence"]
    aln_str = row["Aligned_Sequence"]
    reads = row["#Reads"]

    for sub_ref_pos in row["substitution_positions"]:
        left_index = pos + sub_ref_pos
        char_idx = ref_positions.index(sub_ref_pos)
        yield (chrom, left_index, ref_str[char_idx], aln_str[char_idx], reads)


def build_edit_counts(df_alleles, amplicon_positions):
    """Build a dict mapping (chrom, pos, ref, alt) -> total_reads from allele data.

    Each entry represents one biallelic VCF record. Read counts are
    summed for identical edits across alleles.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        The alleles DataFrame from CRISPResso2 analysis.
    amplicon_positions : dict
        Maps reference name -> (chrom, pos) where pos is the 1-based
        absolute coordinate of reference index 0.

    Returns
    -------
    dict
        Mapping of (chrom, pos, ref, alt) -> total read count.
    """
    edit_counts = defaultdict(int)

    for _, row in df_alleles.iterrows():
        if (
            row.get("n_inserted", 0) == 0
            and row.get("n_deleted", 0) == 0
            and row.get("n_mutated", 0) == 0
        ):
            continue

        ref_name = row["Reference_Name"]
        chrom, pos = amplicon_positions[ref_name]

        for chrom_v, pos_v, ref, alt, reads in _edits_from_deletions(row, chrom, pos):
            edit_counts[(chrom_v, pos_v, ref, alt)] += reads
        for chrom_v, pos_v, ref, alt, reads in _edits_from_insertions(row, chrom, pos):
            edit_counts[(chrom_v, pos_v, ref, alt)] += reads
        for chrom_v, pos_v, ref, alt, reads in _edits_from_substitutions(row, chrom, pos):
            edit_counts[(chrom_v, pos_v, ref, alt)] += reads

    return dict(edit_counts)


def write_vcf_from_edits(edit_counts, num_reads, amplicon_lens, vcf_path):
    """Write biallelic VCF file from edit counts.

    Parameters
    ----------
    edit_counts : dict
        Mapping of (chrom, pos, ref, alt) -> read count.
    num_reads : int
        Total read count (denominator for allele frequencies).
    amplicon_lens : dict
        Maps amplicon/chrom name -> amplicon length (for contig headers).
    vcf_path : str
        Output file path.

    Returns
    -------
    int
        Number of VCF data lines written (excluding header).
    """
    if num_reads is None or num_reads <= 0:
        raise CRISPRessoShared.BadParameterException("Error counting total number of reads.")

    def _key_sort(k):
        chrom, pos, ref, alt = k
        try:
            c_key = (0, int(chrom))
        except (TypeError, ValueError):
            c_key = (1, str(chrom))
        return (c_key, int(pos), len(alt), alt, len(ref), ref)

    counter = 0
    denom = float(num_reads)
    with open(vcf_path, "w") as f:
        f.write('##fileformat=VCFv4.5\n')
        f.write('##source=CRISPResso2\n')
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        for amplicon_name, amplicon_len in amplicon_lens.items():
            f.write(f'##contig=<ID={amplicon_name},length={amplicon_len}>\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key in sorted(edit_counts.keys(), key=_key_sort):
            chrom, pos, ref, alt = key
            af = f"{edit_counts[key] / denom:.3f}"
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tAF={af}\n")
            counter += 1
    return counter


def write_vcf_file(df_alleles, ref_names, ref_lens, args, vcf_path):
    """Orchestrates: parse amplicon coordinates, build edits, and write VCF file.

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
        Solely for testing; the number of VCF data lines written.
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

    num_reads = df_alleles['#Reads'].sum()

    edit_counts = build_edit_counts(df_alleles, amplicon_positions)
    count = write_vcf_from_edits(edit_counts, num_reads, amplicon_lens, vcf_path)

    return count
