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
ctypedef int DTYPE_INT

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


cdef void _diag_cell(size_t k, size_t d, int* LOF, int* OFF,
                     int* dM, int* dI, int* dJ, int* gi, int* sc,
                     int gap_open, int gap_extend,
                     size_t max_i, size_t max_j) noexcept nogil:
    """Compute one interior diagonal cell via the total (i,j)->slot mapping.

    Used for the <=4 edge cells of each diagonal (neighbor windows and the
    gap-rule endpoints j==max_j / i==max_i, which use gap_extend instead of
    gap_open exactly like the row-major last row/col blocks).
    """
    cdef size_t i = <size_t> LOF[d] + k
    cdef size_t j = d - i
    cdef int go = gap_extend if (i == max_i or j == max_j) else gap_open
    cdef size_t dd, i_n
    # left neighbor (i, j-1) on diagonal d-1
    dd = i + j - 1
    i_n = i - <size_t> LOF[dd]
    cdef int a = go + dM[OFF[dd] + i_n] + gi[i]
    cdef int b = gap_extend + dI[OFF[dd] + i_n] + gi[i]
    dI[<size_t> OFF[d] + k] = a if a > b else b
    # up neighbor (i-1, j) on diagonal d-1
    i_n = (i - 1) - <size_t> LOF[dd]
    cdef int c = go + dM[OFF[dd] + i_n] + gi[i - 1]
    cdef int e = gap_extend + dJ[OFF[dd] + i_n]
    dJ[<size_t> OFF[d] + k] = c if c > e else e
    # diagonal neighbor (i-1, j-1) on diagonal d-2
    dd = i + j - 2
    i_n = (i - 1) - <size_t> LOF[dd]
    cdef int m = dM[OFF[dd] + i_n]
    cdef int iv = dI[OFF[dd] + i_n]
    cdef int jv = dJ[OFF[dd] + i_n]
    if iv > m:
        m = iv
    if jv > m:
        m = jv
    dM[<size_t> OFF[d] + k] = m + sc[k]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def global_align(str pystr_seqj, str pystr_seqi, np.ndarray[DTYPE_LONG, ndim=2] matrix,
          np.ndarray[DTYPE_LONG,ndim=1] gap_incentive, int gap_open=-1,
          int gap_extend=-1):
    """
    Global sequence alignment (needleman-wunsch) -- diagonal-major fill.

    Same DP recurrences and edge rules as ``global_align_pointers`` (the
    original row-major implementation, kept as the differential-test
    reference); scores and alignment strings are byte-identical. Two changes
    make the fill ~2-3x faster at long-read sizes:

    * **Diagonal-major storage.** The three score matrices (M/I/J) are stored
      one anti-diagonal at a time in flat buffers with per-diagonal offsets,
      instead of row-major. Cell (i, j) lives at ``BUF[i+j] + i - LO[i+j]``.
      The recurrences' only intra-row serial dependency (I[i,j] needs
      I[i,j-1]) vanishes along an anti-diagonal — every cell's inputs
      (left/up on diagonal d-1, diagonal on d-2) sit at *constant offsets*
      from its own lane, so the per-diagonal inner loop auto-vectorizes
      (NEON) with contiguous loads/stores. Border cells (i==0 or j==0) are
      initialized into their diagonal slots with the same closed-form
      formulas the row-major version writes, making the (i,j) mapping total:
      traceback needs no special cases. Traceback itself walks diagonals in
      decreasing d order, so it reads the buffers near-sequentially.
    * **Pointer-free traceback** (as before): no pointer matrices; the walk
      recomputes each cell's argmax from the score buffers, reproducing the
      original tie-breaking exactly.

    Validated bit-identical to ``global_align_pointers`` by the hypothesis
    differential test (tests/unit_tests/test_CRISPResso2Align.py).
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

    # Near -infinity sentinel for "impossible" border cells.
    cdef int min_score = -1000000000

    # --- setup: int32 matrix/gap-incentive copies, diagonal index tables ---
    mat32n = np.ascontiguousarray(matrix, dtype=np.dtype("i"))
    mat32n = mat32n.reshape(-1)
    gi32n = np.ascontiguousarray(gap_incentive, dtype=np.dtype("i"))
    cdef np.ndarray[DTYPE_INT, ndim=1] mat32 = mat32n
    cdef np.ndarray[DTYPE_INT, ndim=1] gi32 = gi32n
    cdef int* mat = <int*> mat32.data
    cdef int* gi = <int*> gi32.data
    cdef size_t Wm = <size_t> matrix.shape[1]

    cdef size_t D = max_i + max_j          # last diagonal index
    cdef np.ndarray[DTYPE_INT, ndim=1] offn = np.empty(D + 1, dtype=np.dtype("i"))
    cdef np.ndarray[DTYPE_INT, ndim=1] lofn = np.empty(D + 1, dtype=np.dtype("i"))
    cdef int* OFF = <int*> offn.data
    cdef int* LOF = <int*> lofn.data

    # lo_full(d) = max(0, d - max_j); hi_full(d) = min(max_i, d)
    cdef size_t cells = (max_i + 1) * (max_j + 1)
    cdef np.ndarray[DTYPE_INT, ndim=1] dMn = np.empty(cells, dtype=np.dtype("i"))
    cdef np.ndarray[DTYPE_INT, ndim=1] dIn = np.empty(cells, dtype=np.dtype("i"))
    cdef np.ndarray[DTYPE_INT, ndim=1] dJn = np.empty(cells, dtype=np.dtype("i"))
    cdef int* dM = <int*> dMn.data
    cdef int* dI = <int*> dIn.data
    cdef int* dJ = <int*> dJn.data

    cdef size_t d, k, n, o
    o = 0
    for d in range(D + 1):
        OFF[d] = <int> o
        k = d - max_j
        if k > d:  # never (max_j >= 0); keep C-friendly
            k = 0
        if <long> k < 0:
            k = 0
        LOF[d] = <int> k
        # length = min(max_i, d) - k + 1
        n = d if d < max_i else max_i
        n = n - k + 1
        o += n

    # scratch: substitution scores per diagonal (both indices vary per lane)
    cdef size_t cap = (max_i if max_i < max_j else max_j) + 2
    cdef np.ndarray[DTYPE_INT, ndim=1] scn = np.empty(cap, dtype=np.dtype("i"))
    cdef int* sc = <int*> scn.data

    # --- border init (same closed forms as the row-major version) ---
    # row 0: M=min (M[0,0]=0), I=gap_extend*j+gi[0] (I[0,0]=min), J=min
    # col 0: M=min, J=gap_extend*i+gi[0], I=min
    dM[0] = 0
    dI[0] = min_score
    dJ[0] = min_score
    for d in range(1, D + 1):
        o = <size_t> OFF[d]
        if d <= max_j:  # cell (0, d)
            dM[o] = min_score
            dI[o] = gap_extend * <int> d + gi[0]
            dJ[o] = min_score
        if d <= max_i:  # cell (d, 0): last slot of diagonal d
            k = <size_t> ((d if d < max_i else max_i) - LOF[d])
            dM[o + k] = min_score
            dI[o + k] = min_score
            dJ[o + k] = gap_extend * <int> d + gi[0]

    # --- fill, one diagonal at a time ---
    cdef int gi0 = gi[0]
    cdef int a, b, c, e, m, iv, jv
    cdef int i_lo, i_hi, k0, k1, kv0, kv1
    cdef size_t ke0, ks1
    cdef size_t offL, offD, dd
    cdef int *M1p, *I1p, *J1p, *M0p, *I0p, *J0p, *M2p, *I2p, *J2p
    cdef int* gii
    cdef int* gim
    for d in range(2, D + 1):
        # interior cells: i in [max(1, d-max_j), min(max_i, d-1)]
        i_lo = <int> (d - max_j)
        if i_lo < 1:
            i_lo = 1
        i_hi = <int> (d - 1)
        if i_hi > <long> max_i:
            i_hi = <int> max_i
        if i_lo > i_hi:
            continue
        k0 = i_lo - LOF[d]          # interior start index within the diagonal
        k1 = i_hi - LOF[d]
        n = <size_t> (k1 - k0 + 1)
        o = <size_t> OFF[d]
        M2p = dM + o; I2p = dI + o; J2p = dJ + o
        # score gather (per-lane matrix lookup; scalar)
        for k in range(<size_t> k0, <size_t> k1 + 1):
            i = <size_t> (LOF[d] + <int> k)
            sc[k] = mat[<size_t> seqi[i - 1] * Wm + <size_t> seqj[d - i - 1]]
        # scalar edge cells: k in [k0, k0+1] and [k1-1, k1] (neighbor windows),
        # which subsume the <=2 gap-rule endpoint cells (j==max_j at k0 when
        # d > max_j; i==max_i at k1 when d > max_i)
        ke0 = <size_t> (k0 + 1 if k0 + 1 < k1 else k1)
        for k in range(<size_t> k0, ke0 + 1):
            _diag_cell(k, d, LOF, OFF, dM, dI, dJ, gi, sc,
                       gap_open, gap_extend, max_i, max_j)
        ks1 = <size_t> (k1 - 1 if k1 - 1 > k0 + 1 else k0 + 2)
        for k in range(ks1, <size_t> k1 + 1):
            _diag_cell(k, d, LOF, OFF, dM, dI, dJ, gi, sc,
                       gap_open, gap_extend, max_i, max_j)
        # vector middle: constant-offset neighbor loads from diagonals d-1/d-2
        kv0 = k0 + 2
        kv1 = k1 - 2
        if kv1 >= kv0:
            offL = <size_t> (LOF[d] - LOF[d - 1])          # in {0, 1}
            offD = <size_t> (LOF[d] - LOF[d - 2] - 1)      # diag on d-2
            M1p = dM + <size_t> OFF[d - 1]; I1p = dI + <size_t> OFF[d - 1]; J1p = dJ + <size_t> OFF[d - 1]
            M0p = dM + <size_t> OFF[d - 2]; I0p = dI + <size_t> OFF[d - 2]; J0p = dJ + <size_t> OFF[d - 2]
            # k is the ABSOLUTE index within the diagonal (i = LOF[d] + k),
            # so gi[i] = gii[k] with gii based at LOF[d] — not at the
            # interior start i_lo (whose k0 offset differs for d <= max_j).
            gii = gi + <size_t> LOF[d]
            gim = gii - 1
            for k in range(<size_t> kv0, <size_t> kv1 + 1):
                a = gap_open + M1p[k + offL] + gii[k]
                b = gap_extend + I1p[k + offL] + gii[k]
                I2p[k] = a if a > b else b
                c = gap_open + M1p[k + offL - 1] + gim[k]
                e = gap_extend + J1p[k + offL - 1]
                J2p[k] = c if c > e else e
                m = M0p[k + offD]
                iv = I0p[k + offD]
                jv = J0p[k + offD]
                if iv > m:
                    m = iv
                if jv > m:
                    m = jv
                M2p[k] = m + sc[k]

    # --- traceback: identical walk to global_align_pointers, reading scores
    # through the total diagonal mapping X(i,j) = dX[OFF[i+j] + i - LOF[i+j]].
    seqlen = max_i + max_j
    cdef char* tmp_align_j = get_c_string_with_length(seqlen)
    cdef char* tmp_align_i = get_c_string_with_length(seqlen)

    cdef int matchCount = 0
    i = max_i
    j = max_j
    ci = seqi[i - 1] if i > 0 else seqi[0]
    cj = seqj[j - 1] if j > 0 else seqj[0]
    cdef int currMatrix
    currMatrix = MARRAY
    dd = i + j
    m = dM[OFF[dd] + i - LOF[dd]]
    iv = dI[OFF[dd] + i - LOF[dd]]
    jv = dJ[OFF[dd] + i - LOF[dd]]
    if m > jv:
        if m > iv:
            currMatrix = MARRAY
        else:
            currMatrix = IARRAY
    else:
        if jv > iv:
            currMatrix = JARRAY
        else:
            currMatrix = IARRAY

    # Traceback recompute locals. gap_pen is the per-cell gap_open vs gap_extend
    # used during the fill: gap_open only in the strictly-inner region
    # (i < max_i and j < max_j); the last row/col use gap_extend. Matching this
    # exactly is what makes the recomputed argmax reproduce the pointer path.
    cdef int mm, ii2, jj2, gap_pen, jFromMVal, jExtendVal, iFromMVal, iExtendVal
    cdef size_t i_n, j_n
    while i > 0 or j > 0:
        if currMatrix == MARRAY: # came from a match/mismatch: diagonal move
            # which state did the M-cell at (i,j) transition from? argmax of the
            # three predecessor scores at (i-1,j-1); ties resolved toward IARRAY
            # to match the forward fill's nested-if tie-breaking.
            dd = (i - 1) + (j - 1)
            i_n = i - 1 - LOF[dd]
            mm = dM[OFF[dd] + i_n]
            ii2 = dI[OFF[dd] + i_n]
            jj2 = dJ[OFF[dd] + i_n]
            if mm > jj2:
                currMatrix = MARRAY if mm > ii2 else IARRAY
            else:
                currMatrix = JARRAY if jj2 > ii2 else IARRAY
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
            dd = (i - 1) + j
            i_n = i - 1 - LOF[dd]
            jFromMVal = gap_pen + dM[OFF[dd] + i_n] + gi[i - 1]
            jExtendVal = gap_extend + dJ[OFF[dd] + i_n]
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
            dd = i + (j - 1)
            i_n = i - LOF[dd]
            iFromMVal = gap_pen + dM[OFF[dd] + i_n] + gi[i]
            iExtendVal = gap_extend + dI[OFF[dd] + i_n] + gi[i]
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


