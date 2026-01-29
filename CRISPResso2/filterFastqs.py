"""CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
"""

import gzip
import argparse
import io
import os
import datetime
import numpy
import sys


def main():
    """Main entry point for the filterFastqs script."""
    parser = argparse.ArgumentParser(description='Filter fastqs for quality')
    parser.add_argument('--fastq_r1', required=True)
    parser.add_argument('--fastq_r2')
    parser.add_argument('--min_bp_qual_in_read', type=int)
    parser.add_argument('--min_av_read_qual', type=int)
    parser.add_argument('--min_bp_qual_or_N', type=int)
    parser.add_argument('--fastq_r1_out')
    parser.add_argument('--fastq_r2_out')

    args = parser.parse_args()

    filterFastqs(
        args.fastq_r1,
        args.fastq_r2,
        args.fastq_r1_out,
        args.fastq_r2_out,
        args.min_bp_qual_in_read,
        args.min_av_read_qual,
        args.min_bp_qual_or_N,
    )


def filterFastqs(
    fastq_r1=None,
    fastq_r2=None,
    fastq_r1_out=None,
    fastq_r2_out=None,
    min_bp_qual_in_read=None,
    min_av_read_qual=None,
    min_bp_qual_or_N=None,
    debug=False,
):
    """Filter FASTQ files based on quality thresholds for minimum base pair quality, average read quality, and low quality to N conversion.

    Parameters
    ----------
    fastq_r1 : str, optional
        Path to input R1 FASTQ file.
    fastq_r2 : str, optional
        Path to input R2 FASTQ file.
    fastq_r1_out : str, optional
        Path to output R1 FASTQ file.
    fastq_r2_out : str, optional
        Path to output R2 FASTQ file.
    min_bp_qual_in_read : int, optional
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int, optional
        Minimum average read quality threshold.
    min_bp_qual_or_N : int, optional
        Quality threshold below which bases are converted to N.
    debug : bool, optional
        Enable debug output.

    Returns
    -------
    None

    """
    if (debug):
        print('--fastq_r1:' + str(fastq_r1))
        print('--fastq_r2:' + str(fastq_r2))
        print('--min_bp_qual_in_read:' + str(min_bp_qual_in_read))
        print('--min_av_read_qual:' + str(min_av_read_qual))
        print('--min_bp_qual_or_N:' + str(min_bp_qual_or_N))
        print('--fastq_r1_out:' + str(fastq_r1_out))
        print('--fastq_r2_out:' + str(fastq_r2_out))

    startTime = datetime.datetime.now()

    if not os.path.exists(fastq_r1):
        raise Exception("fastq_r1 file '" + fastq_r1 + "' does not exist.")

    if fastq_r2 is not None and not os.path.exists(fastq_r2):
        raise Exception("fastq_r2 file '" + fastq_r2 + "' does not exist.")

    # CREATION OF FILEHANDLES##
    if fastq_r1.endswith('.gz'):
        f1_in = io.BufferedReader(gzip.open(fastq_r1, 'rb'))
        f1_out_filename = fastq_r1.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'
        if fastq_r1_out:
            f1_out_filename = fastq_r1_out

    else:
        f1_in = open(fastq_r1, 'rb')
        f1_out_filename = fastq_r1.replace('.fastq', '') + '_filtered.fastq'
        if fastq_r1_out:
            f1_out_filename = fastq_r1_out

    if f1_out_filename.endswith('.gz'):
        f1_out = gzip.open(f1_out_filename, 'wt')
    else:
        f1_out = open(f1_out_filename, 'w', encoding='utf-8')

    if fastq_r2:
        if fastq_r2.endswith('.gz'):
            f2_in = io.BufferedReader(gzip.open(fastq_r2, 'rb'))
            f2_out_filename = fastq_r2.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'
            if fastq_r2_out:
                f2_out_filename = fastq_r2_out
        else:
            f2_in = open(fastq_r2, 'rb')
            f2_out_filename = fastq_r2.replace('.fastq', '') + '_filtered.fastq'
            if fastq_r2_out:
                f2_out_filename = fastq_r2_out

        if f2_out_filename.endswith('.gz'):
            f2_out = gzip.open(f2_out_filename, 'wt')
        else:
            f2_out = open(f2_out_filename, 'w', encoding='utf-8')
    # END CREATION OF FILEHANDLES##

    if not fastq_r2:
        if min_bp_qual_in_read:
            if min_av_read_qual:
                if min_bp_qual_or_N:
                    run_mBP_mRQ_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
                else:
                    run_mBP_mRQ(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
            elif min_bp_qual_or_N:
                run_mBP_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
            else:
                run_mBP(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        elif min_av_read_qual:
            if min_bp_qual_or_N:
                run_mRQ_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
            else:
                run_mRQ(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        elif min_bp_qual_or_N:
            run_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        else:
            sys.exit('Finished -- No modifications requested')
    elif min_bp_qual_in_read:
        if min_av_read_qual:
            if min_bp_qual_or_N:
                run_mBP_mRQ_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
            else:
                run_mBP_mRQ_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        elif min_bp_qual_or_N:
            run_mBP_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        else:
            run_mBP_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
    elif min_av_read_qual:
        if min_bp_qual_or_N:
            run_mRQ_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
        else:
            run_mRQ_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
    elif min_bp_qual_or_N:
        run_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N)
    else:
        sys.exit('Finished -- No modifications requested')

    endTime = datetime.datetime.now()
    timeDiff = endTime - startTime
    print("Completed in %d seconds\n" % timeDiff.total_seconds())


def run_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by converting low quality bases to N.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
        npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
        f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mRQ(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by minimum average read quality.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        if mean >= min_av_read_qual:
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mBP(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by minimum base pair quality in read.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        if min_qual >= min_bp_qual_in_read:
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mBP_mRQ(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by both minimum base pair quality and average read quality.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        if mean >= min_av_read_qual:
            min_qual = numpy.min(npQualLine)
            if min_qual >= min_bp_qual_in_read:
                f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mBP_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by minimum base pair quality and convert low quality bases to N.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        if min_qual >= min_bp_qual_in_read:
            npSeqLine = numpy.frombuffer(seqLine, 'c')
            npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mRQ_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by minimum average read quality and convert low quality bases to N.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        if mean >= min_av_read_qual:
            npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
            npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


def run_mBP_mRQ_mBPN(f1_in, f1_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter single-end FASTQ by minimum base pair quality, average read quality, and convert low quality bases to N.

    Parameters
    ----------
    f1_in : file
        Input FASTQ file handle.
    f1_out : file
        Output FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        if min_qual >= min_bp_qual_in_read:
            mean = numpy.mean(npQualLine)
            if mean >= min_av_read_qual:
                npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
                npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
                f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')


# PAIRED
def run_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by converting low quality bases to N in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip()
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
        npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
        f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        npSeqLine2 = numpy.frombuffer(seqLine2, 'c').copy()
        npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
        f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, npSeqLine2.tostring().decode('utf-8'), plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))

        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mRQ_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by minimum average read quality in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip().decode('utf-8')
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean2 = numpy.mean(npQualLine2)
        if mean >= min_av_read_qual and mean2 > min_av_read_qual:
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
            f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, seqLine2, plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mBP_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by minimum base pair quality in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip().decode('utf-8')
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual2 = numpy.min(npQualLine2)
        if min_qual >= min_bp_qual_in_read and min_qual2 > min_bp_qual_in_read:
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
            f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, seqLine2, plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))

        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mBP_mRQ_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by both minimum base pair quality and average read quality in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold (unused in this function).

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip().decode('utf-8')
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip().decode('utf-8')
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        mean2 = numpy.mean(npQualLine2)
        if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
            min_qual = numpy.min(npQualLine)
            min_qual2 = numpy.min(npQualLine2)
            if min_qual >= min_bp_qual_in_read and min_qual2 >= min_bp_qual_in_read:
                f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, seqLine, plusLine.decode('utf-8'), qualLine.decode('utf-8')))
                f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, seqLine2, plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mBP_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by minimum base pair quality and convert low quality bases to N in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality (unused in this function).
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip()
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        min_qual2 = numpy.min(npQualLine2)
        if min_qual >= min_bp_qual_in_read and min_qual2 >= min_bp_qual_in_read:
            npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
            npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
            npSeqLine2 = numpy.frombuffer(seqLine2, 'c').copy()
            npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
            f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, npSeqLine2.tostring().decode('utf-8'), plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mRQ_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by minimum average read quality and convert low quality bases to N in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality (unused in this function).
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip()
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        mean = numpy.mean(npQualLine)
        mean2 = numpy.mean(npQualLine2)
        if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
            npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
            npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
            f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
            npSeqLine2 = numpy.frombuffer(seqLine2, 'c').copy()
            npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
            f2_out.write("%s\n%s\n%s\n%s\n" % (idLine2, npSeqLine2.tostring().decode('utf-8'), plusLine2.decode('utf-8'), qualLine2.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


def run_mBP_mRQ_mBPN_pair(f1_in, f1_out, f2_in, f2_out, min_bp_qual_in_read, min_av_read_qual, min_bp_qual_or_N):
    """Filter paired-end FASTQ by minimum base pair quality, average read quality, and convert low quality bases to N in both reads.

    Parameters
    ----------
    f1_in : file
        Input R1 FASTQ file handle.
    f1_out : file
        Output R1 FASTQ file handle.
    f2_in : file
        Input R2 FASTQ file handle.
    f2_out : file
        Output R2 FASTQ file handle.
    min_bp_qual_in_read : int
        Minimum base pair quality threshold for a read to be kept.
    min_av_read_qual : int
        Minimum average read quality threshold.
    min_bp_qual_or_N : int
        Quality threshold below which bases are converted to N.

    """
    idLine = f1_in.readline().rstrip().decode('utf-8')
    idLine2 = f2_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        seqLine2 = f2_in.readline().rstrip()
        plusLine2 = f2_in.readline().rstrip()
        qualLine2 = f2_in.readline().rstrip()

        npQualLine = numpy.frombuffer(qualLine, dtype=numpy.uint8) - 33  # assume illumina 1.7
        npQualLine2 = numpy.frombuffer(qualLine2, dtype=numpy.uint8) - 33  # assume illumina 1.7
        min_qual = numpy.min(npQualLine)
        min_qual2 = numpy.min(npQualLine2)
        if min_qual >= min_bp_qual_in_read and min_qual2 >= min_bp_qual_in_read:
            mean = numpy.mean(npQualLine)
            mean2 = numpy.mean(npQualLine2)
            if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
                npSeqLine = numpy.frombuffer(seqLine, 'c').copy()
                npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
                f1_out.write("%s\n%s\n%s\n%s\n" % (idLine, npSeqLine.tostring().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
                npSeqLine2 = numpy.frombuffer(seqLine2, 'c').copy()
                npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
                f2_out.write(
                    "%s\n%s\n%s\n%s\n"
                    % (idLine2, npSeqLine2.tostring().decode('utf-8'),
                       plusLine2.decode('utf-8'), qualLine2.decode('utf-8'))
                )
        idLine = f1_in.readline().rstrip().decode('utf-8')
        idLine2 = f2_in.readline().rstrip().decode('utf-8')


if __name__ == "__main__":
    # execute only if run as a script
    main()
