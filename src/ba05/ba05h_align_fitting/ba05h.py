"""
Rosalind: BA5H
Find a Highest-Scoring Fitting Alignment of Two Strings

Say that we wish to compare the approximately 20,000 amino acid-long NRP synthetase
from Bacillus brevis with the approximately 600 amino acid-long A-domain from Streptomyces roseosporus,
the bacterium that produces the powerful antibiotic Daptomycin.
We hope to find a region within the longer protein sequence v that has high similarity
with all of the shorter sequence w.
Global alignment will not work because it tries to align all of v to all of w;
local alignment will not work because it tries to align substrings of both v and w.
Thus, we have a distinct alignment application called the Fitting Alignment Problem.

"Fitting" w to v requires finding a substring v' of v
that maximizes the global alignment score between v' and w among all substrings of v.

Fitting Alignment Problem
Construct a highest-scoring fitting alignment between two strings.

Given: Two DNA strings v and w, where v has length at most 10000 and w has length at most 1000.

Return:
The maximum score of a fitting alignment of v and w, followed
by a fitting alignment achieving this maximum score.
Use the simple scoring method in which matches count +1
and both the mismatch and indel penalties are equal to 1.
(If multiple fitting alignments achieving the maximum score exist, you may return any one.)

Sample Dataset
GTAGGCTTAAGGTTA
TAGATA

Sample Output
2
TAGGCTTA
TAGA--TA

═════════════════════════════════════════════════

Info.
  * What is the "fitting alignment"?
    ; region within the longer protein sequence v
      that has high similarity with all of the shorter sequence w
    ; "Fitting" w to v requires finding a substring v' of v
      that maximizes the global alignment score between v' and w among all substrings of v.
    ; local alignment of v, global alignment of w

  * global vs. local vs. fitting (mismatch, indel panelties equal to 1)
        Global                 Local                  Fitting
    v:  GTAGGCTTAAGGTTA        GTAGGCTTAAGGTTA        GTAGGCTTAAGGTTA
    w:   TAGATA                 TAGATA                 TAGATA
        -TAG----A---T-A         TAGATA                 TAGA--TA
    =>  v & w                  v' & w'                v' & w    <-- important
        whole  & whole         substr & substr        substr & whole
        global & global        local  & local         local  & global

Plan 1.
  * dynamic matrix filling
    v :GTAGGCTTAAGGTTA
    w: TAGATA
            T       A       G       A       T       A
        0──-1──-1──-1──-2──-1──-3──-1──-4──-1──-5──-1──-6
        │       │       │       │       │       │       │
    G  -1  -1  -1  -1  -1  +1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1──-2──-1──-1──-1──-2──-1──-3──-1──-4  match, j == 0
        │↖     │       │       │       │       │       │
    T  -1  +1  -1  -1  -1  -1  -1  -1  -1  +1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──+1──-1── 0──-1──-1──-1──-2──-1──-1──-1──-2  match
        │       │↖     │       │       │       │       │
    A  -1  -1  -1  +1  -1  -1  -1  +1  -1  -1  -1  +1  -1
        │       │       │       │       │       │       │
        0──-1── 0──-1──+2──-1──+1──-1── 0──-1──-1──-1── 0  match
        │       │       │↖     │       │       │       │
    G  -1  -1  -1  -1  -1  +1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1──+1──-1──+3──-1──+2──-1──+1──-1── 0  mismatch
        │       │       │       │↖     │       │       │
    G  -1  -1  -1  -1  -1  +1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1── 0──-1──+2──-1──+2──-1──+1──-1── 0  s(i-1,j) == s(i-1,j-1)
        │       │       │       │       ↑       │       │  Insertion(w) takes precedence.
    C  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1──-1──-1──+1──-1──+1──-1──+1──-1── 0  s(i-1,j) == s(i-1,j-1)
        │       │       │       │       ↑       │       │  Insertion(w) takes precedence.
    T  -1  +1  -1  -1  -1  -1  -1  -1  -1  +1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──+1──-1── 0──-1── 0──-1── 0──-1──+2──-1──+1  s(i-1,j) > s(i-1,j-1)
        │↖     │       │       │       │↖     │       │  Match takes precedence.
    T  -1  +1  -1  -1  -1  -1  -1  -1  -1  +1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──+1──-1── 0──-1──-1──-1──-1──-1──+1──-1──+1  s(i-1,j) == s(i-1,j-1)
        │       ↑       │       │       │       │↖     │  Match takes precedence.
    A  -1  -1  -1  +1  -1  -1  -1  +1  -1  -1  -1  +1  -1
        │       │       │       │       │       │       │
        0──-1── 0──-1──+2──-1──+1──-1── 0──-1── 0──-1──[+2] max value
        │       │↖     │       │       │       │       │
    A  -1  -1  -1  +1  -1  -1  -1  +1  -1  -1  -1  +1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1──+1──-1──+1──-1──+2──-1──+1──-1──+1  s(i-1,j-1) == s(i-1,j)
        │       │       │↖     │       │       │       │  Match takes precedence.
    G  -1  -1  -1  -1  -1  +1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1── 0──-1──+2──-1──+1──-1──+1──-1── 0
        │       │       │       │↖     │       │       │
    G  -1  -1  -1  -1  -1  +1  -1  -1  -1  -1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──-1──-1──-1──-1──+1──-1──+1──-1── 0──-1── 0  a(i-1,j) == s(i-1,j-1)
        │       │       │       │       ↑       │       │  Insertion(w) takes precedence.
    T  -1  +1  -1  -1  -1  -1  -1  -1  -1  +1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──+1──-1── 0──-1── 0──-1── 0──-1──+2──-1──+1  s(i-1,j) > s(i-1,j-1)
        │       │       │       │       │↖     │       │  Match takes precedence.
    T  -1  +1  -1  -1  -1  -1  -1  -1  -1  +1  -1  -1  -1
        │       │       │       │       │       │       │
        0──-1──+1──-1── 0──-1──-1──-1──-1──-1──+1──-1──+1
        │       │       │       │       │       │↖     │
    A  -1  -1  -1  +1  -1  -1  -1  +1  -1  -1  -1  +1  -1
        │       │       │       │       │       │       │
        0──-1── 0──-1──+2──-1──+1──-1── 0──-1── 0──-1──[+2]
            T       A       G       A       T       A

  (1) matrix initialization
    - create dynamic table and initialize it
    - types of dynamic table:
      a. dictionary: coordinate-score pairs :: {(int,int):int}
      b. 2-dimensional array
    - two strings
      a. first string : rows (vertical),      local,  substring,      all zeros (TODO: Why?)
      b. second string: columns (horizontal), global, all characters, negative values
  (2) matrix filling
    - nonspecific
  (3) trace back (backtracking)
    - priority checking

═════════════════════════════════════════════════

References:
- CS 466. Introduction to Bioinformatics - Lecture 5
  https://courses.grainger.illinois.edu/cs466/sp2020/slides/Lecture5_WeightedSeqeunceAlignment.pdf
- solution by andreyrozumnyi
  https://github.com/andreyrozumnyi/rosalind/blob/master/Chapter%205/BA5H.py


"""
#!/usr/bin/env python
import time


def lcs_score_backtrack(v, w, penalty):
    """
    (a,a,int) -> {(int,int):int}
    returns a dynamic table (local alignment of v, global alignment of w)
    """
    dt = {}   # dynamic table as dictionary "{(int,int):int}"
    score = 0
    # local alignment of v
    for i in range(len(v)+1):
        dt[(i,0)] = 0
    # global alignment of w
    for j in range(len(w)+1):
        dt[(0,j)] = -penalty * j
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
    returns fitting alignment of two strings
    (local alignment of v, global alignment of w)
    """
    # find max value in column[m], m == len(w)
    max_val = -float('inf')
    for tup, val in dtable.items():
        if val > max_val and tup[1] == len(w):
            max_val = val
    print(max_val)  # <-- print max value
    # find sinks in column[m], maybe multiple sinks
    sinks = set()
    for tup, val in dtable.items():
        if max_val == val and tup[1] == len(w):
            sinks.add(tup)
    # trace back (backtracking)
    aligns = set()
    for sink in sinks:
        new = ([],[])       # new v, w pair by fitting alignment
        i, j = [sink[0], sink[1]]
        while j > 0:
            # match
            if v[i-1] == w[j-1] and dtable[(i,j)] == dtable[(i-1, j-1)] + 1:
                new[0].append(v[i-1])
                new[1].append(w[j-1])
                i -= 1
                j -= 1
            # insertion for w
            elif dtable[(i,j)] == dtable[(i-1,j)] - penalty:
                new[0].append(v[i-1])
                new[1].append("-")
                i -= 1
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
        new = (new[0][::-1], new[1][::-1])
        new = (''.join(new[0]), ''.join(new[1]))
        aligns.add(new)
    return aligns


def main():
    f = open('/home/wsl/rosalind/data/ba05h.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    v = lines[0]
    w = lines[1]
    f.close()

    start_time = time.time()
    dynamic_table = lcs_score_backtrack(v, w, 1)
    answer = lcs_score_output(dynamic_table, v, w, 1)
    i = 1
    for s in answer:
        print("answer", i, ":")
        i += 1
        print(s[0])
        print(s[1])
        print()

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
