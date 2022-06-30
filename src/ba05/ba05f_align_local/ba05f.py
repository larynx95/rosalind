"""
Rosalind: BA5F
Find a Highest-Scoring "Local Alignment" of Two Strings

Local Alignment Problem
Find the highest-scoring local alignment between two strings.

Given: Two amino acid strings.

Return:
The maximum score of a local alignment of the strings,
followed by a local alignment of these strings achieving the maximum score.
Use the PAM250 scoring matrix and indel penalty σ = 5.
(If multiple local alignments achieving the maximum score exist, you may return any one.)

Sample Dataset
MEANLY
PENALTY

Sample Output
15
EANL-Y
ENALTY

* PAM250
PAM250 is one of a family of alignment scoring matrices introduced by Margaret Dayhoff in 1978.
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10
═════════════════════════════════════════════════

Info.
  * global alignment (Needleman-Wunsch)
    F(0,0) = 0
    F(i,j) = max ┌ F(i-1,i-1) + s(x_i, y_j)
                 ├ F(i-1,j) + d
                 └ F(i,j-1) + d

  * local alignment (Smith-Waterman)
    F(0,0) = 0
    F(i,j) = max ┌ F(i-1,j-1) + s(x_i, y_j)
                 ├ F(i-1,j) + d
                 ├ F(i,j-1) + d
                 └ 0

Plan 1.
  * This exercise lacks explanation.
  * Questions:
    - What is the meaning of the phrase, "local alignment of the strings"?
    - What are "global" and "local" alignments?
      global alignment of two strings: entire lengths of two strings
      local alignment of two strings: substrings of two strings
    - How many pairs of substrings can exist? Is it efficient?
      ; so many substring pairs, definitely inefficient
    - What is the "scoring system"?
      ; BLOSUM62, PAM250
    - Which data structure is best for scoring system (BLOSUM62 or PAM250)?
      ; TODO: check this later again

Plan 2.
  * searchig for help in gloogle and youtube
    https://www.youtube.com/watch?v=lu9ScxSejSE
  * Smith-Waterman algorithm
  * example ATGCT, AGCT

    (1) initialization
      m = len(ATGCT) => the number of columns in matrix is m+1
      n = len(AGCT)  => the number of rows    in matrix is n+1

    (2) matrix filling

               A    T    G    C    T
        ----------------------------     match 1
          0    0    0    0    0    0     misamtch -1
        ----------------------------     gap -2
      A   0    1    0    0    0    0     !! all the negative value => zero
        ----------------------------
      G   0    0    0    1    0    0
        ----------------------------
      C   0    0    0    0    2    0
        ----------------------------
      T   0    0    0    0    0    3

    (3) trace back
      highest value == 3, backtrack 3 -> 2 -> 1 TCG => GCT

═════════════════════════════════════════════════

References:
- Smith–Waterman algorithm
  https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
- Local Sequence Alignment & Smith-Waterman || Algorithm and Example
  https://www.youtube.com/watch?v=lu9ScxSejSE
"""
#!/usr/bin/env python
import time


aacids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

dic_blosum62 = {'A': [ 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2],
                'C': [ 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
                'D': [-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3],
                'E': [-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2],
                'F': [-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3],
                'G': [ 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3],
                'H': [-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2],
                'I': [-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1],
                'K': [-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2],
                'L': [-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1],
                'M': [-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1],
                'N': [-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2],
                'P': [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3],
                'Q': [-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1],
                'R': [-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2],
                'S': [ 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2],
                'T': [ 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2],
                'V': [ 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1],
                'W': [-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2],
                'Y': [-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7]}

dic_pam250 = {'A':[ 2,-2, 0, 0,-3, 1,-1,-1,-1,-2,-1, 0, 1, 0,-2, 1, 1, 0,-6,-3],
              'C':[-2,12,-5,-5,-4,-3,-3,-2,-5,-6,-5,-4,-3,-5,-4, 0,-2,-2,-8, 0],
              'D':[ 0,-5, 4, 3,-6, 1, 1,-2, 0,-4,-3, 2,-1, 2,-1, 0, 0,-2,-7,-4],
              'E':[ 0,-5, 3, 4,-5, 0, 1,-2, 0,-3,-2, 1,-1, 2,-1, 0, 0,-2,-7,-4],
              'F':[-3,-4,-6,-5, 9,-5,-2, 1,-5, 2, 0,-3,-5,-5,-4,-3,-3,-1, 0, 7],
              'G':[ 1,-3, 1, 0,-5, 5,-2,-3,-2,-4,-3, 0, 0,-1,-3, 1, 0,-1,-7,-5],
              'H':[-1,-3, 1, 1,-2,-2, 6,-2, 0,-2,-2, 2, 0, 3, 2,-1,-1,-2,-3, 0],
              'I':[-1,-2,-2,-2, 1,-3,-2, 5,-2, 2, 2,-2,-2,-2,-2,-1, 0, 4,-5,-1],
              'K':[-1,-5, 0, 0,-5,-2, 0,-2, 5,-3, 0, 1,-1, 1, 3, 0, 0,-2,-3,-4],
              'L':[-2,-6,-4,-3, 2,-4,-2, 2,-3, 6, 4,-3,-3,-2,-3,-3,-2, 2,-2,-1],
              'M':[-1,-5,-3,-2, 0,-3,-2, 2, 0, 4, 6,-2,-2,-1, 0,-2,-1, 2,-4,-2],
              'N':[ 0,-4, 2, 1,-3, 0, 2,-2, 1,-3,-2, 2, 0, 1, 0, 1, 0,-2,-4,-2],
              'P':[ 1,-3,-1,-1,-5, 0, 0,-2,-1,-3,-2, 0, 6, 0, 0, 1, 0,-1,-6,-5],
              'Q':[ 0,-5, 2, 2,-5,-1, 3,-2, 1,-2,-1, 1, 0, 4, 1,-1,-1,-2,-5,-4],
              'R':[-2,-4,-1,-1,-4,-3, 2,-2, 3,-3, 0, 0, 0, 1, 6, 0,-1,-2, 2,-4],
              'S':[ 1, 0, 0, 0,-3, 1,-1,-1, 0,-3,-2, 1, 1,-1, 0, 2, 1,-1,-2,-3],
              'T':[ 1,-2, 0, 0,-3, 0,-1, 0, 0,-2,-1, 0, 0,-1,-1, 1, 3, 0,-5,-3],
              'V':[ 0,-2,-2,-2,-1,-1,-2, 4,-2, 2, 2,-2,-1,-2,-2,-1, 0, 4,-6,-2],
              'W':[-6,-8,-7,-7, 0,-7,-3,-5,-3,-2,-4,-4,-6,-5, 2,-2,-5,-6,17, 0],
              'Y':[-3, 0,-4,-4, 7,-5, 0,-1,-4,-1,-2,-2,-5,-4,-4,-3,-3,-2, 0,10]}


def graph_to_weighted_edges(graph):
    """
    {a:[a]} -> {(a,a):int}
    returns a weighted edges from a graph (adjacent list)
    """
    aacids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    edges = dict()
    for key, values in graph.items():
        for i in range(len(aacids)):
            edges[(key,aacids[i])] = values[i]
    return edges


edges_pm250 = graph_to_weighted_edges(dic_pam250)


def lcs_score_backtrack(edges, v, w, sigma = 5):
    """
    ({(str,str):int},str,str,int) -> {(int,int):int}
    returns a dynamic table
    """
    s = {}
    score = 0
    for i in range(len(v)+1):
        s[(i,0)] = 0                # s[(i,0)] = -sigma * i
    for j in range(len(w)+1):
        s[(0,j)] = 0                # s[(0,j)] = -sigma * j
    for i in range(len(v)):
        for j in range(len(w)):
            horizontally = 0 if s[(i,j+1)] - sigma < 0 else s[(i,j+1)] - sigma
            vertically   = 0 if s[(i+1,j)] - sigma < 0 else s[(i+1,j)] - sigma
            diagonally   = 0 if s[(i,j)] + edges[(v[i], w[j])] < 0 else s[(i,j)] + edges[(v[i], w[j])]
            s[(i+1,j+1)] = max(horizontally, vertically, diagonally)
    return s


def lcs_score_output(edges, backtrack, v, w, sigma = 5):
    """
    ({(str,str):int},{(int,int):int},str,str,int) -> (str,str)
    returns a local alignment of two strings
    """
    # find max score value
    max_val = -float('inf')
    for tup, score in backtrack.items():
        if score > max_val:
            max_val = score
    print(max_val)
    # find all coordinates with max score - (i,j) tuple
    sinks = set()
    for tup, score in backtrack.items():
        if score == max_val:
            sinks.add(tup)
    # trace back (backtracking)
    aligns = set()
    for sink in sinks:
        new_v = []
        new_w = []
        i, j = [sink[0], sink[1]]
        while backtrack[(i, j)] != 0:
            if backtrack[(i,j)] == backtrack[(i-1,j)] - sigma:
                new_v.append(v[i-1])
                new_w.append('-')
                i -= 1
            elif backtrack[(i,j)] == backtrack[(i,j-1)] - sigma:
                new_v.append('-')
                new_w.append(w[j-1])
                j -= 1
            else:
                new_v.append(v[i-1])
                new_w.append(w[j-1])
                i -= 1
                j -= 1
        new_v, new_w = [new_v[::-1], new_w[::-1]]
        new_v, new_w = [''.join(new_v), ''.join(new_w)]
        aligns.add((new_v, new_w))
    return aligns


def main():
    f = open('/home/wsl/rosalind/data/ba05f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    astr = lines[0]
    bstr = lines[1]
    f.close()

    start_time = time.time()
    backtrack = lcs_score_backtrack(edges_pm250, astr, bstr, 5)
    output = lcs_score_output(edges_pm250, backtrack, astr, bstr, 5)
    for tup in output:
        print(tup[0])
        print(tup[1])
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
