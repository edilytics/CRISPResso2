"""Streaming helpers used by CRISPRessoPooled.

The pooled pipeline still uses samtools as its BAM/FASTA backend, but keeps
record parsing, CIGAR handling, and FASTQ writing in Python.  This avoids the
large embedded shell/awk programs that were previously used for demultiplexing.
"""
from collections import Counter, OrderedDict
from dataclasses import dataclass
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, Iterator, Mapping, Optional, Tuple


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class SamRecord:
    """The SAM fields required to create a FASTQ record."""

    name: str
    flag: int
    reference: str
    position: int
    cigar: str
    sequence: str
    quality: str


def parse_sam_record(line: str) -> Optional[SamRecord]:
    """Parse one non-header SAM line, returning ``None`` for headers."""
    if not line or line.startswith("@"):
        return None
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11:
        raise ValueError(f"Malformed SAM record with {len(fields)} fields: {line!r}")
    return SamRecord(
        name=fields[0],
        flag=int(fields[1]),
        reference=fields[2],
        position=int(fields[3]),
        cigar=fields[5],
        sequence=fields[9],
        quality=fields[10],
    )


def _run_samtools(command, *, stderr=None):
    """Run a samtools command without invoking a shell."""
    return subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=stderr if stderr is not None else subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def iter_sam_records(
    bam_filename: str,
    exclude_flags: str,
    region: Optional[str] = None,
    samtools_command: str = "samtools",
    stderr=None,
) -> Iterator[SamRecord]:
    """Stream alignment records from ``samtools view``.

    ``region`` uses samtools' normal ``chrom:start-end`` syntax.  Errors are
    checked after the stream is consumed so a failed worker cannot silently
    produce an incomplete demultiplexing report.
    """
    command = [samtools_command, "view", "-F", str(exclude_flags), bam_filename]
    if region:
        command.append(region)
    process = _run_samtools(command, stderr=stderr)
    error_output = ""
    return_code = None
    try:
        for line in process.stdout:
            record = parse_sam_record(line)
            if record is not None:
                yield record
    finally:
        process.stdout.close()
        if process.stderr is not None:
            error_output = process.stderr.read()
        return_code = process.wait()
    if return_code:
        detail = error_output.strip()
        raise RuntimeError(f"samtools view failed ({return_code}): {detail}")


def get_bam_references(
    bam_filename: str,
    samtools_command: str = "samtools",
) -> Dict[str, int]:
    """Return reference names and lengths from a BAM header."""
    result = subprocess.run(
        [samtools_command, "view", "-H", bam_filename],
        check=True,
        capture_output=True,
        text=True,
    )
    references = {}
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" in fields and "LN" in fields:
            references[fields["SN"]] = int(fields["LN"])
    return references


def alignment_span(record: SamRecord) -> Tuple[int, int, str]:
    """Return the pooled tool's ``(start, end, strand)`` for a SAM record.

    Coordinates intentionally follow the historical pooled implementation:
    start is one-based and end is exclusive.  Terminal soft clips extend the
    reported span, while insertions and hard clips do not consume reference.
    """
    operations = [(int(length), operation) for length, operation in _CIGAR_RE.findall(record.cigar)]
    start = record.position
    end = record.position
    for length, operation in operations:
        if operation == "S":
            if end == record.position:
                start -= length
            else:
                end += length
        elif operation not in ("I", "H", "P"):
            end += length
    strand = "-" if record.flag & 0x10 else "+"
    return start, end, strand


def fastq_record(record: SamRecord) -> str:
    """Serialize a SAM record as one FASTQ record."""
    return f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n"


class FastqShardWriter:
    """Write many FASTQ shards while bounding open file descriptors."""

    def __init__(self, directory: str, max_open_files: Optional[int] = None):
        self.directory = directory
        self.max_open_files = max_open_files
        self._handles = OrderedDict()
        self._temporary_paths = {}
        self.counts = Counter()
        os.makedirs(directory, exist_ok=True)

    def _temporary_path(self, key) -> str:
        if key not in self._temporary_paths:
            fd, path = tempfile.mkstemp(prefix=".demux_", suffix=".fastq", dir=self.directory)
            os.close(fd)
            self._temporary_paths[key] = path
        return self._temporary_paths[key]

    def _open(self, key):
        if key in self._handles:
            handle = self._handles.pop(key)
            self._handles[key] = handle
            return handle
        if self.max_open_files and len(self._handles) >= self.max_open_files:
            _, handle = self._handles.popitem(last=False)
            handle.close()
        handle = open(self._temporary_path(key), "a", encoding="utf-8")
        self._handles[key] = handle
        return handle

    def write(self, key, record: SamRecord):
        self._open(key).write(fastq_record(record))
        self.counts[key] += 1

    def close(self):
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    def finalize(self, paths: Mapping, minimum_reads: Optional[int] = None, retain_below_threshold: bool = True):
        """Finalize shards and return ``key -> output path or None``."""
        self.close()
        outputs = {}
        for key, count in self.counts.items():
            temporary_path = self._temporary_paths[key]
            final_path = paths[key]
            if minimum_reads is not None and count < minimum_reads:
                if retain_below_threshold:
                    shutil.move(temporary_path, final_path[:-3] if final_path.endswith(".gz") else final_path)
                else:
                    os.remove(temporary_path)
                outputs[key] = None
                continue
            if final_path.endswith(".gz"):
                with open(temporary_path, "rb") as source, gzip.open(final_path, "wb") as destination:
                    shutil.copyfileobj(source, destination)
            else:
                shutil.move(temporary_path, final_path)
            outputs[key] = final_path
        for key, temporary_path in self._temporary_paths.items():
            if os.path.exists(temporary_path):
                os.remove(temporary_path)
        return outputs


def count_fastq_records(fastq_filename: str) -> int:
    """Count FASTQ records without shell utilities."""
    opener = gzip.open if fastq_filename.endswith(".gz") else open
    count = 0
    with opener(fastq_filename, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            record = [header, handle.readline(), handle.readline(), handle.readline()]
            if len(record) != 4 or any(line == "" for line in record[1:]):
                raise ValueError(f"Incomplete FASTQ record in {fastq_filename}")
            count += 1
    return count


def extract_fasta_sequence(reference_filename: str, region: str, samtools_command: str = "samtools") -> str:
    """Extract and concatenate a FASTA region using ``samtools faidx``."""
    result = subprocess.run(
        [samtools_command, "faidx", reference_filename, region],
        check=True,
        capture_output=True,
        text=True,
    )
    return "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">"))


def top_unaligned_sequences(
    bam_filename: str,
    number_of_records: int = 10000,
    number_to_return: int = 10,
    samtools_command: str = "samtools",
) -> list:
    """Return the most common sequences among the first unmapped records."""
    process = _run_samtools([samtools_command, "view", "-f", "4", bam_filename])
    counts = Counter()
    seen = 0
    try:
        for line in process.stdout:
            record = parse_sam_record(line)
            if record is None:
                continue
            counts[record.sequence] += 1
            seen += 1
            if seen >= number_of_records:
                process.terminate()
                break
    finally:
        process.stdout.close()
    if process.stderr is not None:
        process.stderr.close()
    process.wait()
    return [sequence for sequence, _ in counts.most_common(number_to_return)]


def demultiplex_bam_to_fastq(
    bam_filename: str,
    output_paths: Mapping[str, str],
    exclude_flags: str,
    max_open_files: Optional[int] = None,
    samtools_command: str = "samtools",
) -> Dict[str, int]:
    """Demultiplex alignments by reference name into compressed FASTQ files."""
    directory = os.path.dirname(next(iter(output_paths.values()))) if output_paths else "."
    writer = FastqShardWriter(directory, max_open_files=max_open_files)
    for record in iter_sam_records(bam_filename, exclude_flags, samtools_command=samtools_command):
        if record.reference in output_paths:
            writer.write(record.reference, record)
    writer.finalize(output_paths)
    return dict(writer.counts)


def demultiplex_interval(
    bam_filename: str,
    chromosome: str,
    start: int,
    end: int,
    output_directory: str,
    minimum_reads: int,
    exclude_flags: str,
    info_filename: str,
    samtools_command: str = "samtools",
) -> str:
    """Write one restricted genomic interval and its depth manifest."""
    os.makedirs(output_directory, exist_ok=True)
    stem = os.path.join(output_directory, f"REGION_{chromosome}_{start}_{end}.fastq")
    output_path = stem + ".gz"
    writer = FastqShardWriter(output_directory, max_open_files=1)
    key = (chromosome, start, end)
    for record in iter_sam_records(
        bam_filename,
        exclude_flags,
        region=f"{chromosome}:{start}-{end}",
        samtools_command=samtools_command,
    ):
        writer.write(key, record)
    count = writer.counts.get(key, 0)
    writer.finalize({key: output_path}, minimum_reads=minimum_reads, retain_below_threshold=True)
    reported_path = output_path if count >= minimum_reads else "NA"
    with open(info_filename, "w", encoding="utf-8") as handle:
        handle.write(f"{chromosome}\t{start}\t{end}\t{count}\t{reported_path}\n")
    return info_filename


def demultiplex_jobs_chunk(jobs):
    """Run a chunk of pooled demultiplexing jobs in a worker process."""
    results = []
    for job in jobs:
        if job["mode"] == "interval":
            results.append(demultiplex_interval(**job["arguments"]))
        else:
            results.append(demultiplex_genome_chunk(**job["arguments"]))
    return results


def demultiplex_genome_chunk(
    bam_filename: str,
    chromosome: str,
    region_start: Optional[int],
    region_end: Optional[int],
    output_directory: str,
    minimum_reads: int,
    exclude_flags: str,
    info_filename: str,
    samtools_command: str = "samtools",
) -> str:
    """Group a chromosome/chunk by calculated alignment span."""
    os.makedirs(output_directory, exist_ok=True)
    region = None if region_start is None else f"{chromosome}:{region_start}-{region_end}"
    writer = FastqShardWriter(output_directory, max_open_files=32)
    for record in iter_sam_records(
        bam_filename,
        exclude_flags,
        region=region,
        samtools_command=samtools_command,
    ):
        start, end, _ = alignment_span(record)
        writer.write((record.reference, start, end), record)

    paths = {
        key: os.path.join(output_directory, f"REGION_{key[0]}_{key[1]}_{key[2]}.fastq.gz")
        for key in writer.counts
    }
    writer.finalize(paths, minimum_reads=minimum_reads, retain_below_threshold=True)
    with open(info_filename, "w", encoding="utf-8") as handle:
        for (reference, start, end), count in sorted(writer.counts.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])):
            output_path = paths[(reference, start, end)] if count >= minimum_reads else "NA"
            handle.write(f"{reference}\t{start}\t{end}\t{count}\t{output_path}\n")
    return info_filename
