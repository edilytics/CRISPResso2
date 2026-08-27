import gzip
import stat

from CRISPResso2 import CRISPRessoPooledDemux as demux


SAM = """read1\t0\tAMPL_a\t1\t255\t2S4M1I2M1S\t*\t0\t0\tAACCGGT\tIIIIIII
read2\t16\tAMPL_b\t3\t255\t5M\t*\t0\t0\tACGTA\tIIIII
"""


def fake_samtools(tmp_path, output=SAM):
    script = tmp_path / "samtools"
    script.write_text("#!/bin/sh\ncat <<'EOF'\n" + output + "EOF\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return str(script)


def test_alignment_span_handles_soft_clipping_and_strand():
    record = demux.parse_sam_record(SAM.splitlines()[0])
    assert demux.alignment_span(record) == (-1, 8, "+")

    reverse_record = demux.parse_sam_record(SAM.splitlines()[1])
    assert demux.alignment_span(reverse_record) == (3, 8, "-")


def test_fastq_shard_writer_bounds_open_files_and_compresses(tmp_path):
    writer = demux.FastqShardWriter(str(tmp_path), max_open_files=1)
    first = demux.parse_sam_record(SAM.splitlines()[0])
    writer.write("a", first)
    writer.write("b", first)
    paths = {"a": str(tmp_path / "a.fastq.gz"), "b": str(tmp_path / "b.fastq.gz")}
    writer.finalize(paths)

    with gzip.open(paths["a"], "rt") as handle:
        assert handle.read().startswith("@read1\n")
    assert writer.counts == {"a": 1, "b": 1}


def test_demultiplex_bam_to_fastq_uses_python_records(tmp_path):
    samtools = fake_samtools(tmp_path)
    paths = {
        "AMPL_a": str(tmp_path / "AMPL_a.fastq.gz"),
        "AMPL_b": str(tmp_path / "AMPL_b.fastq.gz"),
    }
    counts = demux.demultiplex_bam_to_fastq("reads.bam", paths, "0x900", samtools_command=samtools)
    assert counts == {"AMPL_a": 1, "AMPL_b": 1}
    with gzip.open(paths["AMPL_b"], "rt") as handle:
        assert handle.read() == "@read2\nACGTA\n+\nIIIII\n"


def test_top_unaligned_sequences(tmp_path):
    unmapped = """r1\t4\t*\t0\t0\t*\t*\t0\t0\tAAA\tIII
r2\t4\t*\t0\t0\t*\t*\t0\t0\tCCC\tIII
r3\t4\t*\t0\t0\t*\t*\t0\t0\tAAA\tIII
"""
    samtools = fake_samtools(tmp_path, unmapped)
    assert demux.top_unaligned_sequences("reads.bam", 10, 2, samtools) == ["AAA", "CCC"]


def test_count_fastq_records_and_extract_fasta(tmp_path):
    fastq = tmp_path / "reads.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        handle.write("@r\nAC\n+\nII\n")
    assert demux.count_fastq_records(str(fastq)) == 1

    script = tmp_path / "samtools-faidx"
    script.write_text("#!/bin/sh\nprintf '>chr1:1-4\\nAC\\nGT\\n'")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    assert demux.extract_fasta_sequence("ref.fa", "chr1:1-4", str(script)) == "ACGT"


def test_genome_chunk_restricts_chromosome_wide_jobs(tmp_path):
    args_file = tmp_path / "samtools.args"
    script = tmp_path / "samtools"
    script.write_text(
        f"#!/bin/sh\nprintf '%s\\n' \"$@\" > '{args_file}'\n"
        "cat <<'EOF'\n"
        + SAM
        + "EOF\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)

    demux.demultiplex_genome_chunk(
        "reads.bam", "chr1", None, None, str(tmp_path), 1, "0x900", str(tmp_path / "chunk.info"), str(script)
    )

    assert args_file.read_text().splitlines()[-1] == "chr1"


def test_interval_and_genome_chunk_reports(tmp_path):
    samtools = fake_samtools(tmp_path)
    mapped = tmp_path / "mapped"
    interval_info = mapped / "interval.info"
    demux.demultiplex_interval(
        "reads.bam", "chr1", 1, 10, str(mapped), 1, "0x900", str(interval_info), samtools_command=samtools
    )
    assert "chr1\t1\t10\t2\t" in interval_info.read_text()

    chunk_info = mapped / "chunk.info"
    demux.demultiplex_genome_chunk(
        "reads.bam", "chr1", None, None, str(mapped), 1, "0x900", str(chunk_info), samtools_command=samtools
    )
    lines = chunk_info.read_text().splitlines()
    assert any(line.startswith("AMPL_a\t-1\t8\t1\t") for line in lines)
    assert any(line.startswith("AMPL_b\t3\t8\t1\t") for line in lines)
