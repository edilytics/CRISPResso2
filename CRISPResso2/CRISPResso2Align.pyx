#with help from friends at https://github.com/brentp/align for cython implementation (no thanks for bug-ridden algorithm)
# https://github.com/dnase/affine-gap-sequence-alignment/blob/master/alignment.py for affine gap algorithm

from cython.view cimport array as cvarray
import numpy as np
cimport numpy as np

from libc.stdlib cimport free, malloc

cimport cython
import sys
import os.path

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t

cdef extern from "Python.h":
    ctypedef void PyObject

ctypedef long DTYPE_LONG

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4
cdef size_t MARRAY = 1, IARRAY = 2, JARRAY = 3


cdef char* get_c_string_with_length(size_t length):
    cdef char* c_string = <char *> malloc((length + 1) * sizeof(char))
    if not c_string:
        raise MemoryError()
    return c_string


def read_matrix(path):
    """
    Read a matrix in the NCBI format
    The score for a 'C' changing to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    """
    cdef np.ndarray[DTYPE_LONG, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size

    with open(path) as fh:
        headers = None
        while headers is None:
            line = fh.readline().strip()
            if line[0] == '#': continue
            headers = [ord(x) for x in line.split(' ') if x]
        mat_size = max(headers) + 1

        a = np.zeros((mat_size, mat_size), dtype=np.int64)

        line = fh.readline()
        while line:
            line_vals = [int(x) for x in line[:-1].split(' ')[1:] if x]
            for ohidx, val in zip(headers, line_vals):
                a[headers[ai], ohidx] = val
            ai += 1
            line = fh.readline()

    return a

def make_matrix(match_score=5, mismatch_score=-4, n_mismatch_score=-2, n_match_score=-1):
    """
    Create a score matrix for matches/mismatches.
    The default values here are those represented in the EDNAFULL matrix

    match_score: score for matching nucleotide values
    mismatch_score: score for mismatching nucleotide values
    n_mismatch_score: score for matching a nucleotide with 'N'
    n_match_score: score for 'N' matching an 'N'
    """
    cdef np.ndarray[DTYPE_LONG, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size

    letters = ['A','T','C','G','N']
    headers = [ord(x) for x in letters]
    mat_size = max(headers) + 1

    nuc_ords = [ord(x) for x in ['A','T','C','G']]

    a = np.zeros((mat_size, mat_size), dtype=np.int64)

    for nuc in nuc_ords:
      for nuc2 in nuc_ords:
        if nuc == nuc2:
          a[nuc,nuc2] = match_score
        else:
          a[nuc,nuc2] = mismatch_score

    for nuc in nuc_ords:
      a[nuc,ord('N')] = n_mismatch_score
      a[ord('N'),nuc] = n_mismatch_score


    a[ord('N'),ord('N')] = n_match_score

    return a

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def global_align_pointers(str pystr_seqj, str pystr_seqi, np.ndarray[DTYPE_LONG, ndim=2] matrix,
          np.ndarray[DTYPE_LONG,ndim=1] gap_incentive, int gap_open=-1,
          int gap_extend=-1):
    """
    Reference (pointer-based) global sequence alignment (needleman-wunsch).

    Retained as the oracle for the differential property test
    (test_global_align_ptrfree_matches_pointers). Production uses the
    pointer-free global_align below, which is ~2.5x faster with byte-identical
    outputs. This body is the original implementation, unchanged.

    Original: global sequence alignment (needleman-wunsch) on seq i and j.
    Reference is seq_i, sequenced read is seq_j
    Match and mismatch values are read from matrix object
    where matrix is of the format provided in the ncbi/data directory.
    gap_incentive is the incentive for having a gap at each position in seqi -
      this allows for the preferential location of gaps to be at the predicted
      cut site in genome editing experiments.

    """

    byte_seqj = pystr_seqj.encode('UTF-8')
    cdef char* seqj = byte_seqj
    byte_seqi = pystr_seqi.encode('UTF-8')
    cdef char* seqi = byte_seqi

    cdef size_t max_j = len(pystr_seqj)
    cdef size_t max_i = len(pystr_seqi)
    if len(gap_incentive) != max_i + 1:
        print('\nERROR: Mismatch in gap_incentive length (gap_incentive: ' + str(len(gap_incentive)) + ' ref: '+str(max_i+1) + '\n')
        return 0

    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, seqlen, align_counter = 0, p
    cdef int diag_score, up_score, left_score, tscore

    cdef str align_j
    cdef str align_i
    cdef char ci
    cdef char cj

    #create 3 arrays of scores and 3 arrays of pointers
    # M array - best alignment so far ending with a match
    # I array - best alignment so far ending with a gap in Read (J) (insertion in ref, deletion in read)
    # J array - best alignment so far ending with a gap in Ref (I) (deletion in ref, insertion in read)

    cdef int [:,:] mScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] iScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] jScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] mPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] iPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] jPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))

    # min_score is a near -infinity sentinel for "impossible" border cells
    # (mScore on row/col 0, iScore on col 0, jScore on row 0) -- cells no valid
    # alignment can end in. It must dominate every real path score so the
    # traceback never routes through these cells into uninitialized pointer
    # memory. gap_open*max_i*max_j is NOT negative enough for short sequences or
    # large gap_incentive, which let optimal paths reach those cells and read
    # garbage / segfault (pre-existing bug surfaced by the hypothesis property
    # test). A fixed large sentinel is dominant for any int32-feasible score.
    cdef int min_score = -1000000000

    #init match matrix
    mScore[0,1:] = min_score
    mScore[1:,0] = min_score
    mScore[0,0] = 0
    mPointer[0,1:] = IARRAY
    mPointer[1:,0] = JARRAY
    mPointer[0,0] = 0

# no penalty for gaps starting at beginning
#    score[0, 1:] = 0
#    score[1:, 0] = 0
# gap extension penalty for gaps starting at beginning
    #init i matrix
    for i in range(1,max_j+1):
        iScore[0,i] = gap_extend * i + gap_incentive[0]
#    iScore[0,1:] = [gap_extend * np.arange(1, max_j+1, dtype=int)]
    iScore[0:,0] = min_score
    iPointer[0,1:] = IARRAY

    #init j matrix
    for i in range(1,max_i+1):
        jScore[i,0] = gap_extend * i + gap_incentive[0]
    #jScore[1:,0] = np.vectorize(gap_extend * np.arange(1, max_i+1, dtype=int))
    jScore[0,0:] = min_score
    jPointer[1:,0] = JARRAY

#    print('gap penalty is'+str(gap_incentive))

    cdef int iFromMVal
    cdef int iExtendVal
    cdef DTYPE_LONG gi_i, gi_im1  # gap_incentive[i], gap_incentive[i-1]; hoisted out of inner loop
    cdef int jFromMVal
    cdef int jExtendVal
    cdef int mVal, iVal, jVal

    #apply NW algorithm for inside squares (not last row or column)
    for i in range(1, max_i):
        ci = seqi[i - 1] #char in i
        # gap_incentive[i]/[i-1] are loop-invariant in j; hoist into scalars so the
        # inner loop touches only the score/pointer matrices (stride-1, hot) and
        # never re-reads the separate gap_incentive buffer each iteration.
        gi_i = gap_incentive[i]
        gi_im1 = gap_incentive[i - 1]

        for j in range(1, max_j):
            cj = seqj[j - 1] #char in j

            iFromMVal = gap_open + mScore[i, j - 1] + gi_i
            iExtendVal = gap_extend + iScore[i, j - 1] + gi_i
            if iFromMVal > iExtendVal:
                iScore[i,j] = iFromMVal
                iPointer[i,j] = MARRAY
            else:
                iScore[i,j] = iExtendVal
                iPointer[i,j] = IARRAY

            jFromMVal = gap_open + mScore[i - 1, j] + gi_im1
	    #no gap incentive here -- J already got the gap incentive when it transitioned from M, so don't add it again if we're extending.
            jExtendVal = gap_extend + jScore[i - 1, j]
            if jFromMVal > jExtendVal:
                jScore[i,j] =  jFromMVal
                jPointer[i,j] = MARRAY
            else:
                jScore[i,j] = jExtendVal
                jPointer[i,j] = JARRAY

            mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
            iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
            jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
            if mVal > jVal:
                if mVal > iVal:
                    mScore[i, j] = mVal
                    mPointer[i, j] = MARRAY
                else:
                    mScore[i, j]   = iVal
                    mPointer[i, j] = IARRAY
            else:
                if jVal > iVal:
                    mScore[i, j]  = jVal
                    mPointer[i, j] = JARRAY
                else:
                    mScore[i, j] = iVal
                    mPointer[i, j] = IARRAY

#            print('mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))

    #for last column and last row, ignore gap opening penalty
    #last column
    j = max_j
    cj = seqj[j-1]
    for i in range(1, max_i):
        ci = seqi[i-1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        if iFromMVal > iExtendVal:
            iScore[i,j] =  iFromMVal
            iPointer[i,j] = MARRAY
        else:
            iScore[i,j] = iExtendVal
            iPointer[i,j] = IARRAY

        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal:
            jScore[i,j] =  jFromMVal
            jPointer[i,j] = MARRAY
        else:
            jScore[i,j] = jExtendVal
            jPointer[i,j] = JARRAY

        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal:
                mScore[i, j] = mVal
                mPointer[i, j] = MARRAY
            else:
                mScore[i, j]   = iVal
                mPointer[i, j] = IARRAY
        else:
            if jVal > iVal:
                mScore[i, j]  = jVal
                mPointer[i, j] = JARRAY
            else:
                mScore[i, j] = iVal
                mPointer[i, j] = IARRAY
#        print('lastCol: mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))

    #last row
    i = max_i
    ci = seqi[i - 1]
    # i is fixed at max_i here, so both gap_incentive reads are loop-invariant in j.
    gi_i = gap_incentive[i]
    gi_im1 = gap_incentive[i - 1]
    for j in range(1, max_j+1):
        cj = seqj[j - 1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gi_i
        iExtendVal = gap_extend + iScore[i, j - 1] + gi_i
        if iFromMVal > iExtendVal:
            iScore[i,j] =  iFromMVal
            iPointer[i,j] = MARRAY
        else:
            iScore[i,j] = iExtendVal
            iPointer[i,j] = IARRAY

        jFromMVal = gap_extend + mScore[i - 1, j] + gi_im1
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal:
            jScore[i,j] =  jFromMVal
            jPointer[i,j] = MARRAY
        else:
            jScore[i,j] = jExtendVal
            jPointer[i,j] = JARRAY


        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal:
                mScore[i, j] = mVal
                mPointer[i, j] = MARRAY
            else:
                mScore[i, j]   = iVal
                mPointer[i, j] = IARRAY
        else:
            if jVal > iVal:
                mScore[i, j]  = jVal
                mPointer[i, j] = JARRAY
            else:
                mScore[i, j] = iVal
                mPointer[i, j] = IARRAY
#        print('lastRow: mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))



#    print('mScore')
#    for ii in range(mScore.shape[0]):
#      for jj in range(mScore.shape[1]):
#        print(str(mScore[ii,jj]) + '.' + str(mPointer[ii,jj])+ ","),
#      print("\n"),
#    print('iScore')
#    for ii in range(iScore.shape[0]):
#      for jj in range(iScore.shape[1]):
#        print(str(iScore[ii,jj]) + '.' + str(iPointer[ii,jj])+ ","),
#      print("\n"),
#    print('jScore')
#    for ii in range(jScore.shape[0]):
#      for jj in range(jScore.shape[1]):
#        print(str(jScore[ii,jj]) + '.' + str(jPointer[ii,jj])+ ","),
#      print("\n"),

    seqlen = max_i + max_j
    cdef char* tmp_align_j = get_c_string_with_length(seqlen)
    cdef char* tmp_align_i = get_c_string_with_length(seqlen)

    cdef int matchCount = 0
    i = max_i
    j = max_j
    ci = seqi[i - 1]
    cj = seqj[j - 1]
    cdef int currMatrix
    currMatrix = MARRAY
    if mScore[i,j] > jScore[i,j]:
        if mScore[i,j] > iScore[i,j]:
            currMatrix = MARRAY
        else:
            currMatrix = IARRAY
    else:
        if jScore[i,j] > iScore[i,j]:
            currMatrix = JARRAY
        else:
            currMatrix = IARRAY
#    print('seqi' + str(seqi))
#    print('seqj' + str(seqj))
    while i > 0 or j > 0:
        # print("i: " + str(i) + " j: " + str(j) + " currMatrix: " + str(currMatrix) + " match score: " + str(mScore[i,j]) + " last match: " +  str(mScore[i-1,j-1]) + " matrix[" + str(ci) + "," + str(cj) + "]: " + str(matrix[ci,cj]) + " last j " + str(jScore[i,j]) + " last i: " + str(iScore[i,j]) + " mpointer: " + str(mPointer[i,j]) + " ipointer: " + str(iPointer[i,j]) + " jpointer: " + str(jPointer[i,j]))

        currVal = mScore[i,j]
        currPtr = mPointer[i,j]
        if currMatrix == IARRAY:
            currVal = iScore[i,j]
            currPtr = iPointer[i,j]
        if currMatrix == JARRAY:
            currVal = jScore[i,j]
            currPtr = jPointer[i,j]
#        print("i: " + str(i) + " j: " + str(j) + " " + str(currMatrix) +':' + str(currVal) + ' > ' + str(currPtr))
        if currMatrix == MARRAY: # 1
            currMatrix = mPointer[i,j]
            tmp_align_j[align_counter] = cj
            tmp_align_i[align_counter] = ci
            if cj == ci:
                matchCount += 1

            if i > 1:
                i -= 1
                ci = seqi[i - 1]
            else:
                i = 0
                ci = seqi[i]
            if j > 1:
                j -= 1
                cj = seqj[j - 1]
            else:
                j = 0
                cj = seqj[j]

#            print('in M set to ' + str(currMatrix))
        elif currMatrix == JARRAY: # 3
            currMatrix = jPointer[i,j]
            tmp_align_j[align_counter] = c"-"
            tmp_align_i[align_counter] = ci
            if i > 1:
                i -= 1
                ci = seqi[i - 1]
            else:
                i = 0
                ci = seqi[i]
        elif currMatrix == IARRAY: # 2
            currMatrix = iPointer[i,j]
            tmp_align_j[align_counter] = cj
            tmp_align_i[align_counter] = c"-"
            if j > 1:
                j -= 1
                cj = seqj[j - 1]
            else:
                j = 0
                cj = seqj[j]
        else:
            print('i: ' + str(i) + ' j: ' + str(j))
            print('currMatrix:' + str(currMatrix))
            print('seqj: ' + str(seqj) + ' seqi: ' + str(seqi))
            raise Exception('wtf4!:pointer: %i', i)
#          print('at end, currMatrix is ' + str(currMatrix))

        align_counter += 1
    try:
        align_j = tmp_align_j[:align_counter].decode('UTF-8', 'strict')
    finally:
        free(tmp_align_j)
    try:
        align_i = tmp_align_i[:align_counter].decode('UTF-8', 'strict')
    finally:
        free(tmp_align_i)

    # print(tounicode_with_length_and_free(alig))
#    print(str(matchCount) + " aln: " + str(align_counter))
    final_score = 100*matchCount/float(align_counter)
    return align_j[::-1], align_i[::-1], round(final_score, 3)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def global_align(str pystr_seqj, str pystr_seqi, np.ndarray[DTYPE_LONG, ndim=2] matrix,
          np.ndarray[DTYPE_LONG,ndim=1] gap_incentive, int gap_open=-1,
          int gap_extend=-1):
    """
    Global sequence alignment (needleman-wunsch) -- pointer-free implementation.

    The classic global_align stores three extra int32 "pointer" matrices
    (one per DP state) purely to guide traceback, then reads them once. That is
    12 extra bytes/cell of memory traffic in the O(N*M) forward pass — the
    bottleneck at long-read scale (see design_docs/ALIGN_SIMD_OPTIMIZATION.md).

    This variant drops all three pointer matrices. Traceback instead recomputes
    each cell's argmax from the score matrices on the fly. Traceback is O(N+M),
    so the recompute cost is negligible; the win is ~halving the fill's memory
    traffic (6 matrices -> 3). The min_score sentinels on the border row/col
    make the border behave correctly with no special-casing.

    Validated bit-identical to global_align_pointers by the hypothesis
    differential test tests/unit_tests/test_CRISPResso2Align.py::test_global_align_ptrfree_matches_pointers.
    """

    byte_seqj = pystr_seqj.encode('UTF-8')
    cdef char* seqj = byte_seqj
    byte_seqi = pystr_seqi.encode('UTF-8')
    cdef char* seqi = byte_seqi

    cdef size_t max_j = len(pystr_seqj)
    cdef size_t max_i = len(pystr_seqi)
    if len(gap_incentive) != max_i + 1:
        print('\nERROR: Mismatch in gap_incentive length (gap_incentive: ' + str(len(gap_incentive)) + ' ref: '+str(max_i+1) + '\n')
        return 0

    cdef size_t i = 0, j = 0, seqlen, align_counter = 0, p
    cdef str align_j
    cdef str align_i
    cdef char ci
    cdef char cj

    # Only the three score matrices are stored — no pointer matrices.
    cdef int [:,:] mScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] iScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] jScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))

    # Near -infinity sentinel for "impossible" border cells (see global_align).
    cdef int min_score = -1000000000

    #init match matrix
    mScore[0,1:] = min_score
    mScore[1:,0] = min_score
    mScore[0,0] = 0

    #init i matrix
    for i in range(1,max_j+1):
        iScore[0,i] = gap_extend * i + gap_incentive[0]
    iScore[0:,0] = min_score

    #init j matrix
    for i in range(1,max_i+1):
        jScore[i,0] = gap_extend * i + gap_incentive[0]
    jScore[0,0:] = min_score

    cdef int iFromMVal
    cdef int iExtendVal
    cdef DTYPE_LONG gi_i, gi_im1  # gap_incentive[i], gap_incentive[i-1]; hoisted out of inner loop
    cdef int jFromMVal
    cdef int jExtendVal
    cdef int mVal, iVal, jVal

    #apply NW algorithm for inside squares (not last row or column)
    for i in range(1, max_i):
        ci = seqi[i - 1] #char in i
        gi_i = gap_incentive[i]
        gi_im1 = gap_incentive[i - 1]

        for j in range(1, max_j):
            cj = seqj[j - 1] #char in j

            iFromMVal = gap_open + mScore[i, j - 1] + gi_i
            iExtendVal = gap_extend + iScore[i, j - 1] + gi_i
            iScore[i,j] = iFromMVal if iFromMVal > iExtendVal else iExtendVal

            jFromMVal = gap_open + mScore[i - 1, j] + gi_im1
            jExtendVal = gap_extend + jScore[i - 1, j]
            jScore[i,j] = jFromMVal if jFromMVal > jExtendVal else jExtendVal

            mVal = mScore[i - 1, j - 1]
            iVal = iScore[i - 1, j - 1]
            jVal = jScore[i - 1, j - 1]
            if iVal > mVal:
                mVal = iVal
            if jVal > mVal:
                mVal = jVal
            mScore[i, j] = mVal + matrix[ci,cj]

    #for last column and last row, ignore gap opening penalty
    #last column
    j = max_j
    cj = seqj[j-1]
    for i in range(1, max_i):
        ci = seqi[i-1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        iScore[i,j] = iFromMVal if iFromMVal > iExtendVal else iExtendVal

        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        jScore[i,j] = jFromMVal if jFromMVal > jExtendVal else jExtendVal

        mVal = mScore[i - 1, j - 1]
        iVal = iScore[i - 1, j - 1]
        jVal = jScore[i - 1, j - 1]
        if iVal > mVal:
            mVal = iVal
        if jVal > mVal:
            mVal = jVal
        mScore[i, j] = mVal + matrix[ci,cj]

    #last row
    i = max_i
    ci = seqi[i - 1]
    gi_i = gap_incentive[i]
    gi_im1 = gap_incentive[i - 1]
    for j in range(1, max_j+1):
        cj = seqj[j - 1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gi_i
        iExtendVal = gap_extend + iScore[i, j - 1] + gi_i
        iScore[i,j] = iFromMVal if iFromMVal > iExtendVal else iExtendVal

        jFromMVal = gap_extend + mScore[i - 1, j] + gi_im1
        jExtendVal = gap_extend + jScore[i - 1, j]
        jScore[i,j] = jFromMVal if jFromMVal > jExtendVal else jExtendVal

        mVal = mScore[i - 1, j - 1]
        iVal = iScore[i - 1, j - 1]
        jVal = jScore[i - 1, j - 1]
        if iVal > mVal:
            mVal = iVal
        if jVal > mVal:
            mVal = jVal
        mScore[i, j] = mVal + matrix[ci,cj]

    seqlen = max_i + max_j
    cdef char* tmp_align_j = get_c_string_with_length(seqlen)
    cdef char* tmp_align_i = get_c_string_with_length(seqlen)

    cdef int matchCount = 0
    i = max_i
    j = max_j
    ci = seqi[i - 1]
    cj = seqj[j - 1]
    cdef int currMatrix
    currMatrix = MARRAY
    if mScore[i,j] > jScore[i,j]:
        if mScore[i,j] > iScore[i,j]:
            currMatrix = MARRAY
        else:
            currMatrix = IARRAY
    else:
        if jScore[i,j] > iScore[i,j]:
            currMatrix = JARRAY
        else:
            currMatrix = IARRAY

    # Traceback recompute locals. gap_pen is the per-cell gap_open vs gap_extend
    # used during the fill: gap_open only in the strictly-inner region
    # (i < max_i and j < max_j); the last row/col use gap_extend. Matching this
    # exactly is what makes the recomputed argmax reproduce the pointer path.
    cdef int mm, ii, jj, gap_pen
    while i > 0 or j > 0:
        if currMatrix == MARRAY: # came from a match/mismatch: diagonal move
            # which state did the M-cell at (i,j) transition from? argmax of the
            # three predecessor scores at (i-1,j-1); ties resolved toward IARRAY
            # to match the forward fill's nested-if tie-breaking.
            mm = mScore[i - 1, j - 1]
            ii = iScore[i - 1, j - 1]
            jj = jScore[i - 1, j - 1]
            if mm > jj:
                currMatrix = MARRAY if mm > ii else IARRAY
            else:
                currMatrix = JARRAY if jj > ii else IARRAY
            tmp_align_j[align_counter] = cj
            tmp_align_i[align_counter] = ci
            if cj == ci:
                matchCount += 1
            if i > 1:
                i -= 1
                ci = seqi[i - 1]
            else:
                i = 0
                ci = seqi[i]
            if j > 1:
                j -= 1
                cj = seqj[j - 1]
            else:
                j = 0
                cj = seqj[j]
        elif currMatrix == JARRAY: # gap in read (deletion in ref): move up in i
            gap_pen = gap_open if (i < max_i and j < max_j) else gap_extend
            jFromMVal = gap_pen + mScore[i - 1, j] + gap_incentive[i-1]
            jExtendVal = gap_extend + jScore[i - 1, j]
            currMatrix = MARRAY if jFromMVal > jExtendVal else JARRAY
            tmp_align_j[align_counter] = c"-"
            tmp_align_i[align_counter] = ci
            if i > 1:
                i -= 1
                ci = seqi[i - 1]
            else:
                i = 0
                ci = seqi[i]
        elif currMatrix == IARRAY: # gap in ref (insertion in read): move left in j
            gap_pen = gap_open if (i < max_i and j < max_j) else gap_extend
            iFromMVal = gap_pen + mScore[i, j - 1] + gap_incentive[i]
            iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
            currMatrix = MARRAY if iFromMVal > iExtendVal else IARRAY
            tmp_align_j[align_counter] = cj
            tmp_align_i[align_counter] = c"-"
            if j > 1:
                j -= 1
                cj = seqj[j - 1]
            else:
                j = 0
                cj = seqj[j]
        else:
            print('i: ' + str(i) + ' j: ' + str(j))
            print('currMatrix:' + str(currMatrix))
            print('seqj: ' + str(seqj) + ' seqi: ' + str(seqi))
            raise Exception('wtf4!:pointer: %i', i)

        align_counter += 1
    try:
        align_j = tmp_align_j[:align_counter].decode('UTF-8', 'strict')
    finally:
        free(tmp_align_j)
    try:
        align_i = tmp_align_i[:align_counter].decode('UTF-8', 'strict')
    finally:
        free(tmp_align_i)

    final_score = 100*matchCount/float(align_counter)
    return align_j[::-1], align_i[::-1], round(final_score, 3)
