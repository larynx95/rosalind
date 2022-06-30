"""
Rosalind: BA5E
Find a Highest-Scoring "Global Alignment" of Two Strings

Global Alignment Problem
Find the highest-scoring alignment between two strings using a scoring matrix.

Given:
Two amino acid strings.

Return:
The maximum alignment score of these strings followed by an alignment achieving this maximum score.
Use the BLOSUM62 scoring matrix and indel penalty σ (sigma) = 5.
(If multiple alignments achieving the maximum score exist, you may return any one.)

Sample Dataset
PLEASANTLY
MEANLY

Sample Output
8
PLEASANTLY
-MEA--N-LY

* Additional information: BLOSUM62
BLOSUM62, given in detail below,
is a commonly used scoring matrix in alignment problems considering protein strings.
It assigns different alignment scores to substituted amino acids
depending on the particular acids substituted;
scores were computed based on the relative frequency of one amino acid to be substituted for another
in a collection of known alignments.

The scoring matrix is "symmetric",
which in this case means that the score assigned to substituting symbol x for symbol y
is the same as that of substituting y for x.

  |  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
--|------------------------------------------------------------
A |  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C |  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D | -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E | -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F | -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G |  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H | -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I | -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K | -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L | -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M | -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N | -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P | -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q | -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R | -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S |  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T |  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V |  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W | -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y | -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7

═════════════════════════════════════════════════


Plan 1.
  * using algorithm in BA5C
    rectangular grid, (a) right, (b) down, (c) diagonal

  * sample dataset: score 8
    PLEASANTLY   (-: insertion)
    -MEA--N-LY   (-: deletion)     ↘: match or mismatch, →: insertion, ↓: deletion

          M   E   A   N   L   Y
        ┌───┬───┬───┬───┬───┬───┐
      P ↓   │   │   │   │   │   │  deletion  -5
        ├───┼───┼───┼───┼───┼───┤
      L │↘ │   │   │   │   │   │  mismatch   2
        ├───┼───┼───┼───┼───┼───┤
      E │   │↘ │   │   │   │   │  match      5
        ├───┼───┼───┼───┼───┼───┤
      A │   │   │↘ │   │   │   │  match      4
        ├───┼───┼───┼───┼───┼───┤
      S │   │   │   ↓   │   │   │  deletion  -5
        ├───┼───┼───┼───┼───┼───┤
      A │   │   │   ↓   │   │   │  deletion  -5
        ├───┼───┼───┼───┼───┼───┤
      N │   │   │   │↘ │   │   │  match      6
        ├───┼───┼───┼───┼───┼───┤
      T │   │   │   │   ↓   │   │  deletion  -5
        ├───┼───┼───┼───┼───┼───┤
      L │   │   │   │   │↘ │   │  match      4
        ├───┼───┼───┼───┼───┼───┤
      Y │   │   │   │   │   │↘ │  match      7  total 8
        └───┴───┴───┴───┴───┴───┘

Plan 2.
  * "Needleman-Wunsch algorithm"
    https://www.youtube.com/watch?v=ipp-pNRIp4g
    https://www.youtube.com/watch?v=um8h3P216Fk
    In this algorithm, data structure is matrix not adjacency list.

  * example   ATGCT, AGCT

    (1) initialization
      m = len(ATGCT) => the number of columns in matrix is m+1
      n = len(AGCT)  => the number of rows    in matrix is n+1

    (2) matrix filling

               A    T    G    C    T
        ----------------------------     match 1
          0   -2   -4   -6   -8  -10     misamtch -1
        ----------------------------     gap -2
      A  -2    1   -1   -3   -5   -7
        ----------------------------
      G  -4   -1    0    0   -2   -4
        ----------------------------
      C  -6   -3   -2   -1    1   -1
        ----------------------------
      T  -8   -5   -2   -3   -1    2

    (3) backtrack

═════════════════════════════════════════════════

References:
- Global Sequence Alignment & Needleman-Wunsch || Algorithm and Example
  https://www.youtube.com/watch?v=ipp-pNRIp4g
"""
#!/usr/bin/env python
import time


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


edges_blosum62 = graph_to_weighted_edges(dic_blosum62)


def lcs_score_backtrack(edges, v, w, sigma):
    """
    ({(str,str):int},str,str,int) -> {(int,int):int}
    returns a dynamic table
    >>> lcs_score_backtrack(edges_blosum62, 'PLEASANTLY','MEANLY')
        {(0,0):0, (1,0):-5, ... , (10,6):8}
    """
    s = {}
    score = 0
    for i in range(len(v)+1):
        s[(i,0)] = -sigma * i
    for j in range(len(w)+1):
        s[(0,j)] = -sigma * j
    for i in range(len(v)):
        for j in range(len(w)):
            horizontally = s[(i,j+1)] - sigma
            vertically   = s[(i+1,j)] - sigma
            diagonally   = s[(i,j)] + edges[(v[i], w[j])]
            s[(i+1,j+1)] = max(horizontally, vertically, diagonally)
    print(s[(len(v),len(w))])
    return s


def lcs_score_output(edges, backtrack, v, w, sigma):
    """
    ({(str,str):int},{(int,int):int},str,str,int) -> (str,str)
    returns a global alignment of two strings
    """
    i = len(v)
    j = len(w)
    new_v = []
    new_w = []
    while i != 0 or j != 0:  # 'or' !!
        if backtrack[(i,j)] == backtrack[(i-1,j)] - sigma:
            new_v.append(v[i-1])
            new_w.append('-')
            i -= 1
        elif backtrack[(i,j)] == backtrack[(i,j-1)] - sigma:
            new_v.append('-')
            new_w.append(w[j-1])
            j -= 1
        else:
        #elif backtrack[(i,j)] == backtrack[(i-1,j-1)] + edges[(v[i-1], w[j-1])]:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
            i -= 1
            j -= 1
    if i != 0:
        new_w += ['-'] * i
    if j != 0:
        new_v += ['-'] * j
    new_v, new_w = [new_v[::-1], new_w[::-1]]
    return ''.join(new_v), ''.join(new_w)


def main():
    f = open('/home/wsl/rosalind/data/ba05e.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    str_a = lines[0]
    str_b = lines[1]
    f.close()

    start_time = time.time()
    backtrack = lcs_score_backtrack(edges_blosum62, str_a, str_b, 5)
    output = lcs_score_output(edges_blosum62, backtrack, str_a, str_b, 5)
    print(output[0])
    print(output[1])
    backtrack = lcs_score_backtrack(edges_blosum62, 'PRTEINS', 'PRTWPSEIN', 5)
    output = lcs_score_output(edges_blosum62, backtrack, 'PRTEINS', 'PRTWPSEIN', 5)
    print(output[0])
    print(output[1])
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
