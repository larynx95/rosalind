"""
Rosalind: BA5I
Find a Highest-Scoring Overlap Alignment of Two Strings

When we assembled genomes, we discussed how to use overlapping reads to assemble a genome,
a problem that was complicated by errors in reads.
We would like to find overlaps between error-prone reads as well.

An overlap alignment of strings v = v1 ... vn and w = w1 ... wm
is a global alignment of a "suffix of v" with a "prefix of w".
An optimal overlap alignment of strings v and w maximizes the global alignment score
between an "i-suffix of v" and a "j-prefix of w"
(i.e., between vi ... vn and w1 ... wj) among all i and j.

Overlap Alignment Problem
Construct a highest-scoring overlap alignment between two strings.

Given: Two protein strings v and w, each of length at most 1000.

Return: The score of an optimal overlap alignment of v and w,
followed by an alignment of a suffix v' of v and a prefix w' of w achieving this maximum score.
Use an alignment score in which matches count +1 and both the mismatch and indel penalties are 2.
(If multiple overlap alignments achieving the maximum score exist, you may return any one.)

Sample Dataset
PAWHEAE
HEAGAWGHEE

Sample Output
1
HEAE
HEAG

═════════════════════════════════════════════════

Info.
  * What is the "overlap alignment"?
    ; "global" alignment of a suffix of v with a prefix of w
       ------                 ------             ------  <- substrings
       ^ Is this really global?
    ; "local" alignment for v and w
      "global" alignment for a suffix of v with a prefix of w

Plan 1.
  * steps:
    (1) local alignment of two strings
    (2) max score, find a suffix of v, prefix of w
        extract left lower corner of matrix rectangle

  * matrix filling
    - left lower corner of matrix rectangle
    - I can see a little bit after drawing the matrix.
      However, I am not sure why local alignment should be applied to two strings.
    - TODO: What exactly do the words local" and "global" mean?
            And what does that have to do with initializing the matrix?

         H       E       A       G       A       W       G       H       E       E  <= prefix of this
     0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0──-2── 0
     │       │       │       │       │       │       │       │       │       │       │
  P -2  -2  -2  -2  -2  -2  -2  -2  -2      -2      -2      -2      -2      -2      -2
     │       │       │       │       │       │       │       │       │       │       │
     0──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2──-2
     │       │       │       │       │       │       │       │       │       │       │
  A -2  -2  -2  -2  -2  +1  -2  -2  -2      -2      -2      -2      -2      -2      -2
     │       │       │       │       │       │       │       │       │       │       │
     0──-2──-2──-2──-4──-2──-1──-2──-3──-2──-1──-2──-3──-2──-4──-2──-4──-2──-4──-2──-4
     │       │       │       │       │       │       │       │       │       │       │
  W -2  -2  -2  -2  -2  -2  -2  -2  -2      -2      -2      -2      -2      -2      -2
     │       │       │       │       │       │       │       │       │       │       │
     0══-2══-2══-2══-4══-2══-3══-2══-3──-2──-3──-2── 0──-2──-2──-2──-4──-2──-6──-2──-6
     ║↖     ║       ║       ║       ║       │       │       │       │       │       │
  H -2  +1  -2  -2  -2  -2  -2  -2  -2      -2      -2      -2      -2      -2      -2
     ║       ║       ║       ║       ║       │       │       │       │       │       │
     0══-2══+1══-2══-1══-2══-3══-2══-5──-2──-5──-2──-2──-2──-2──-2──-1──-2──-3──-2──-5
     ║       ║↖     ║       ║       ║       │       │       │       │       │       │
  E -2  -2  -2  +1  -2  -2  -2  -2  -2      -2      -2      -2      -2      -2      -2
     ║       ║       ║       ║       ║       │       │       │       │       │       │
     0══-2══-1══-2══+2══-2══ 0══-2══-2──-2──-4──-2──-4──-2──-4──-2──-3──-2── 0──-2──-2
     ║       ║       ║↖     ║       ║       │       │       │       │       │       │
  A -2  -2  -2  -2  -2  +1  -2  -2  -2      -2      -2      -2      -2      -2      -2
     ║       ║       ║       ║       ║       │       │       │       │       │       │
     0══-2══-2══-2══ 0══-2══+3══-2══+1──-2──-1──-2──-3──-2──-5──-2──-5──-2──-2──-2──-2
     ║       ║       ║       ║↖     ║       │       │       │       │       │       │
  E -2  -2  -2  +1  -2  -2  -2  -2  -2      -2      -2      -2      -2      -2      -2
     ║       ║       ║       ║       ║       │       │       │       │       │       │
     0══-2══-2══-2══-1══-2═[+1]═-2─[+1]─-2──-1──-2──-3──-2──-5──-2──-7──-2──-4──-2──-1
                             ^       ^
                             max value

═════════════════════════════════════════════════

References:
"""

#!/usr/bin/env python
import time


def lcs_score_backtrack(v, w, penalty=2):
    """
    (a,a,int) -> {(int,int):int}
    returns a dynamic table (local alignment)
    """
    dt = {}   # dynamic table as dictionary "{(int,int):int}"
    score = 0
    # local alignment of v
    for i in range(len(v)+1):
        dt[(i,0)] = 0
    # local alignment of w
    for j in range(len(w)+1):
        dt[(0,j)] = 0
    # create dynamic table as dictionary, not matrix (2d array)
    for i in range(len(v)):
        for j in range(len(w)):
            horizontally = dt[(i,j+1)] - penalty
            vertically   = dt[(i+1,j)] - penalty
            diagonally   = dt[(i,j)] + 1 if v[i] == w[j] else dt[(i,j)] - penalty
            dt[(i+1,j+1)] = max(horizontally, vertically, diagonally)
    return dt


def lcs_score_output(dtable, v, w, penalty):
    """
    ({(int,int):int},a,a,int) -> (a,a)
    returns overlap alignment of two strings
    """
    # find max value in column[m], m == len(w)
    max_val = -float('inf')
    for tup, val in dtable.items():
        if val > max_val and tup[0] == len(v):
            max_val = val
    print(max_val)  # <-- print max value
    # find sinks in column[m], maybe multiple sinks
    sinks = set()
    for tup, val in dtable.items():
        if max_val == val and tup[0] == len(v):
            sinks.add(tup)
    # trace back (backtracking)
    aligns = set()
    for sink in sinks:
        new = ([],[])       # new v, w pair by fitting alignment
        i, j = [sink[0], sink[1]]
        while j > 0 and i > 0:
            # match
            if v[i-1] == w[j-1] and dtable[(i,j)] == dtable[(i-1, j-1)] + 1:
                new[0].append(v[i-1])
                new[1].append(w[j-1])
                i -= 1
                j -= 1
            # deletion for w
            elif dtable[(i,j)] == dtable[(i,j-1)] - penalty:
                new[0].append('-')
                new[1].append(w[j-1])
                j -= 1
            # mismatch
            elif dtable[(i,j)] == dtable[(i-1,j-1)] - penalty:
                new[0].append(v[i-1])
                new[1].append(w[j-1])
                i -= 1
                j -= 1
            # insertion for w
            elif dtable[(i,j)] == dtable[(i-1,j)] - penalty:
                new[0].append(v[i-1])
                new[1].append("-")
                i -= 1
        new = (new[0][::-1], new[1][::-1])
        new = (''.join(new[0]), ''.join(new[1]))
        aligns.add(new)
    return aligns


def main():
    f = open('/home/wsl/rosalind/data/ba05i.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    v = lines[0]
    w = lines[1]
    f.close()

    start_time = time.time()
    backtrack = lcs_score_backtrack(v, w, 2)
    answers = lcs_score_output(backtrack, v, w, 2)
    for ans in answers:
        print(ans[0])
        print(ans[1])
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
