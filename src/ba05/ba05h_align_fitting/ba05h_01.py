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

Given:
Two DNA strings v and w, where v has length at most 10000 and w has length at most 1000.

Return:
The maximum score of a fitting alignment of v and w,
followed by a fitting alignment achieving this maximum score.
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

References:
-

"""

import sys


def fitting_alignment(s1, s2):
    n, m = len(s1), len(s2)
    d = [[0] * (m + 1) for _ in range(n + 1)]
    # local for s1, global for s2
    for j in range(1, m + 1):
        d[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            q = 1 if s1[i-1] == s2[j-1] else -1
            d[i][j] = max(d[i-1][j] - 1, d[i][j-1] - 1, d[i-1][j-1] + q)
    print(d)
    # searching for the best score
    index = 0
    s = d[0][m]
    for i in range(1, n):
        if d[i][m] > s:
            s = d[i][m]
            index = i
    i, j = index, m
    r1, r2 = '', ''
    # backtracking
    while j > 0 and i > 0:
        if d[i][j] == d[i - 1][j - 1] + 1 and s1[i - 1] == s2[j - 1]:
            r1 += s1[i - 1]
            r2 += s2[j - 1]
            i -= 1
            j -= 1
        elif d[i][j] == d[i][j - 1] - 1:
            r1 += '-'
            r2 += s2[j - 1]
            j -= 1
        elif d[i][j] == d[i - 1][j - 1] - 1:
            r1 += s1[i - 1]
            r2 += s2[j - 1]
            i -= 1
            j -= 1
        elif d[i][j] == d[i - 1][j] - 1:
            r1 += s1[i - 1]
            r2 += '-'
            i -= 1
    return d[index][m], r1[::-1], r2[::-1]


def main():
    f = open('/home/wsl/rosalind/data/ba05h.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    v = lines[0]
    w = lines[1]
    f.close()

    for o in fitting_alignment(v, w):
        print(o)


if __name__ == '__main__':
    main()
