"""
Rosalind: BA5C
Find a Longest Common Subsequence of Two Strings

Longest Common Subsequence Problem
Given: Two strings.

Return: A longest common subsequence of these strings.

Sample Dataset
AACCTTGG
ACACTGTGA

Sample Output
AACTGG

═════════════════════════════════════════════════



Plan 1.
  * using algorithms in textbook
  ╔════════════════════════════════════════════════════════════════════════════════╗
  ║ LONGESTPATH(Graph, source, sink)                                               ║
  ║     for each node b in Graph                                                   ║
  ║     s_b <- - infinite                                                          ║
  ║     s_source <- 0                                                              ║
  ║     topologically order Graph   # <-- already topological in this exercise     ║
  ║     for each node b in Graph (following the topological order)                 ║
  ║     s_b <- max<all predecessors a of node b>{s_a + weight of edge from a to b} ║
  ║     return s_sink                                                              ║
  ╚════════════════════════════════════════════════════════════════════════════════╝

  ╔═════════════════════════════════════════════════════════════════╗
  ║ LCSBACKTRACK(v, w)                                              ║
  ║     for i <- 0 to |v|                                           ║
  ║         s(i,o) <- 0                                             ║
  ║     for j <- 0 to |w|                                           ║
  ║         s(0,j) <- 0                                             ║
  ║     for i <- 1 to |v|                                           ║
  ║         for j <- 1 to |w|                                       ║
  ║                           ┌ s(i-1,j)                            ║
  ║             s(i,j) <- max ┼ s(i,j-1)                            ║
  ║                           └ s(i-1,j-1) + 1, if v_i == v_j       ║
  ║             if s(i,j) == s(i-1.j)                               ║
  ║                 Backtrack(i,j) <- "↓"                           ║
  ║             else if s(i,j) == s(i,j-1)                          ║
  ║                 Backtrack(i,j) <- "→"                           ║
  ║             else if s(i,j) == s(i-1,j-1) + 1 and v_i == w_j     ║
  ║                 Backtrack(i,j) <- "↘"                          ║
  ║     return Backtrack                                            ║
  ╚═════════════════════════════════════════════════════════════════╝

  ╔════════════════════════════════════════════╗
  ║ OUTPUTLCS(Backtrack, v, i, j)              ║
  ║     if i = 0 or j = 0                      ║
  ║         return                             ║
  ║     if Backtrack(i,j) = ↓                  ║
  ║         OUTPUTLCS(Backtrack, v, i-1, j)    ║
  ║     else if Backtrack(i,j) = →             ║
  ║         OUTPUTLCS(Backtrack, v, i, j-1)    ║
  ║     else                                   ║
  ║         OUTPUTLCS(Backtrack, v, i-1, j-1)  ║
  ║         output v_i                         ║
  ╚════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════════════╗
  ║ TOPOLOGICALORDERING(Graph)                                           ║
  ║     List <- empty list                                               ║
  ║     Candidates <- set of all nodes in Graph with no incoming edges   ║
  ║     while Candidates is non-empty                                    ║
  ║         select an arbitrary node a from Candidates                   ║
  ║         add a to the end of List and remove it from Candidates       ║
  ║         for each outgoing edge from a to another node b              ║
  ║             remove edge (a, b) from Graph                            ║
  ║             if b has no incoming edges                               ║
  ║                 add b to Candidates                                  ║
  ║     if Graph has edges that have not been removed                    ║
  ║         return "the input graph is not a DAG"                        ║
  ║     else                                                             ║
  ║         return List                                                  ║
  ╚══════════════════════════════════════════════════════════════════════╝

  * This is very similar to Manhattan tourist problem!
  * sample dataset (match/mismatch: ↘/↘, insertion: →, deletion: ↓)
    These two strings were aligned by the solution of Manhattan tourist problem.
    Without it, I can't even align these two strings.
    A T - G T T A T A        (-: insertion)           weight 0
    A T C G T - C - C        (-: deletion)            weight 0
    ↘↘→ ↘↘↓ ↘↓ ↘       (↘: match or mismatch)  weight (score) 1 if two characters are the same

    Growing alignment   Remaingin symbols  Score
    --------------------------------------------
                        A T G T T A T A
                          A T C G T C C

    A                     T G T T A T A    +1
    A                       T C G T C C

    A T                     G T T A T A    +1
    A T                       C G T C C

    A T -                   G T T A T A
    A T C                       G T C C

    A T - G                   T T A T A    +1
    A T C G                       T C C

    A T - G T                   T A T A    +1
    A T C G T                       C C

    A T - G T T                   A T A
    A T C G T -                     C C

    A T - G T T A                   T A
    A T C G T - C                     C

    A T - G T T A T                   A
    A T C G T - C -                   C

    A T - G T T A T A
    A T C G T - C - C
    * *   * *

          A   T   C   G   T   C   C
        ┌───┬───┬───┬───┬───┬───┬───┐
      A │↘ │   │   │   │   │   │   │   match      +1
        ├───┼───┼───┼───┼───┼───┼───┤
      T │   │↘ │   │   │   │   │   │   match      +1
        ├───┼───┼─→─┼───┼───┼───┼───┤   insertion
      G │   │   │   │↘ │   │   │   │   match      +1
        ├───┼───┼───┼───┼───┼───┼───┤
      T │   │   │   │   │↘ │   │   │   match      +1
        ├───┼───┼───┼───┼───┼───┼───┤
      T │   │   │   │   │   ↓   │   │   deletion
        ├───┼───┼───┼───┼───┼───┼───┤
      A │   │   │   │   │   │↘ │   │   mismatch
        ├───┼───┼───┼───┼───┼───┼───┤
      T │   │   │   │   │   │   ↓   │   deletion
        ├───┼───┼───┼───┼───┼───┼───┤
      A │   │   │   │   │   │   │↘ │   mismatch
        └───┴───┴───┴───┴───┴───┴───┘

  * another sample dataset
    a. case 1

    A - A C C T T G - G -
    A C A C - T - G T G A
    ↘→ ↘↘↓ ↘↓ ↘→ ↘→
    A   A C   T   G   G     -->  AACTGG

    Growing alignment      Remaingin symbols      Score
    ---------------------------------------------------
                             A A C C T T G G
                           A C A C T G T G A

    A                          A C C T T G G
    A                        C A C T G T G A

    A -                        A C C T T G G
    A C                        A C T G T G A

    A - A                        C C T T G G
    A C A                        C T G T G A

    A - A C                        C T T G G
    A C A C                        T G T G A

    A - A C C                        T T G G
    A C A C -                      T G T G A

    A - A C C T                        T G G
    A C A C - T                      G T G A

    A - A C C T T                        G G
    A C A C - T -                    G T G A

    A - A C C T T G                        G
    A C A C - T - G                    T G A

    A - A C C T T G -                      G
    A C A C - T - G T                    G A

    A - A C C T T G - G
    A C A C - T - G T G                    A

    A - A C C T T G - G -
    A C A C - T - G T G A
    *   * *   *   *   *

          A   C   A   C   T   G   T   G   A
        ┌───┬───┬───┬───┬───┬───┬───┬───┬───┐
      A │↘ │   │   │   │   │   │   │   │   │  match     +1
        ├───┼─→─┼───┼───┼───┼───┼───┼───┼───┤  insertion
      A │   │   │↘ │   │   │   │   │   │   │  match     +1
        ├───┼───┼───┼───┼───┼───┼───┼───┼───┤
      C │   │   │   │↘ │   │   │   │   │   │  match     +1
        ├───┼───┼───┼───┼───┼───┼───┼───┼───┤
      C │   │   │   │   ↓   │   │   │   │   │  deletion
        ├───┼───┼───┼───┼───┼───┼───┼───┼───┤
      T │   │   │   │   │↘ │   │   │   │   │  match     +1
        ├───┼───┼───┼───┼───┼───┼───┼───┼───┤
      T │   │   │   │   │   ↓   │   │   │   │  deletion
        ├───┼───┼───┼───┼───┼───┼───┼───┼───┤
      G │   │   │   │   │   │↘ │   │   │   │  match     +1
        ├───┼───┼───┼───┼───┼───┼─→─┼───┼───┤  insertion
      G │   │   │   │   │   │   │   │↘ │   │  match     +1
        └───┴───┴───┴───┴───┴───┴───┴───┴─→─┘  insertion

    b. case 2

      A - A C C T - T G G
      A C A C - T G T G A
      ↘→ ↘↘↓ ↘→ ↘↘↘
      A   A C   T   T G     -->  AACTTG

  * rectangular matrix from two strings --> already topologically ordered
    only two arrows (↓, →), no diagonal arrow, as in ManhattanTourist problem

═════════════════════════════════════════════════

References:
"""

import time


#################################################
# Review Manhattan
#################################################


def manhattan_tourist(n, m, down, right):
    """
    (int,int,[[int]],[[int]]) -> int
    algorithm in textbook
    using dynamic programming, directed acyclic graph (DAG)
    >>> manhattan_tourist(4,4,[[1,0,2,4,3],[4,6,5,2,1],[4,4,5,2,1],[5,6,8,5,3]],
                              [[3,2,4,0],[3,2,4,2],[0,7,3,3],[3,3,0,2],[1,3,2,2]])
        34
    TODO: Compare this function with LCS function below.
          The two functions are very similar to each other.
    """
    nodes = [[0] * (m+1) for i in range(n+1)]
    for i in range(1,n+1):
        nodes[i][0] = nodes[i-1][0] + down[i-1][0]
    for j in range(1,m+1):
        nodes[0][j] = nodes[0][j-1] + right[0][j-1]
    for i in range(1,n+1):
        for j in range(1,m+1):
            nodes[i][j] = max(nodes[i-1][j] + down[i-1][j], nodes[i][j-1] + right[i][j-1])
    return nodes[n][m]


#################################################
# backtrack + output = lcs
#################################################


def lcs_backtrack(v,w):
    """
    (str,str) -> set{(int,int)}
    >>> lcs_backtrack("AACT","ACACT")
        {(0,0):0,(1,0):0,(2,0):0,(3,0):0,(4,0):0,
         (0,1):0,(0,2):0,(0,3):0,(0,4):0,(0,5):0,
         (1,1):1,(1,2):1,(1,3):1,(1,4):1,(1,5):1,
         (2,1):1,(2,2):1,(2,3):2,(2,4):2,(2,5):2,
         (3,1):1,(3,2):2,(3,3):2,(3,4):3,(3,5):3,
         (4,1):1,(4,2):2,(4,3):2,(4,4):3,(4,5):4}
    """
    s = {}
    for i in range(len(v)+1):
        s[(i,0)]= 0
    for j in range(len(w)+1):
        s[(0,j)] = 0
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                s[(i+1,j+1)] = s[(i,j)] + 1
            else:
                candidates = [s[(i+1,j)],s[(i,j+1)]]
                s[(i+1,j+1)] = max(candidates)
    return s


def output_lcs(backtrack,v,w):
    """
    ({(int,int)},str,str) -> str
    >>> output_lcs(lcs_backtrack("AACT","ACACT"),"AACT","ACACT")
        AACT
    """
    i = len(v)
    j = len(w)
    result = []
    while i*j !=0:
        if backtrack[(i,j)] == backtrack[(i-1,j)]:
            i -= 1
        elif backtrack[(i,j)] == backtrack[(i,j-1)]:
            j -= 1
        else:
            result.append(v[i-1])
            j -= 1
            i -= 1
    return ''.join(result[::-1])


#################################################
# lcs = backtrack + output
#################################################


def longest_common_subsequence(str1, str2):
    """
    (str,str) -> str
    dynamic programming
    >>> longest_common_subsequence('AACCTTGG', 'ACACTGTGA')
        AACTGG
    """
    # constructing matrix
    m, n = len(str1), len(str2)
    matrix = [[0]*(n+1) for i in range(m+1)]
    for i in range(m+1):
        for j in range(n+1):
            if i==0 or j==0:
                matrix[i][j] = 0
            elif str1[i-1] == str2[j-1]:
                matrix[i][j] = matrix[i-1][j-1]+1
            else:
                matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1])
    # backtracking
    res = ""
    i, j = len(str1), len(str2)
    while (i>0 and j>0):
        if matrix[i][j] == matrix[i-1][j]:
            i -= 1
        elif matrix[i][j] == matrix[i][j-1]:
            j -= 1
        else:
            res = str1[i-1]+res
            i -= 1
            j -= 1
    print(matrix[m][n])
    return res


def longest_common_subsequence2(v, w):
    v = '-' + v
    w = '-' + w

    S = [[0 for _ in range(len(w))] for _ in range(len(v))]
    Backtrack = [[None for _ in range(len(w))] for _ in range(len(v))]

    for i in range(1, len(v)):
        for j in range(1, len(w)):
            tmp = S[i - 1][j - 1] + (1 if v[i] == w[j] else 0)
            S[i][j] = max(S[i - 1][j], S[i][j - 1], tmp)

            if S[i][j] == S[i - 1][j]:
                Backtrack[i][j] = "up"
            elif S[i][j] == S[i][j - 1]:
                Backtrack[i][j] = "left"
            else:
                Backtrack[i][j] = "diag"
    LCS = ""
    while i > 0 and j > 0:
        if Backtrack[i][j] == "diag":
            LCS = v[i] + LCS
            i -= 1
            j -= 1
        elif Backtrack[i][j] == "left":
            j -= 1
        else:
            i -= 1
    return LCS


def main():
    try:
        with open('/home/wsl/rosalind/data/ba05c.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            str1 = lines[0]
            str2 = lines[1]
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    print(longest_common_subsequence('PLEASANTLY','MEANLY'))
    print(longest_common_subsequence2('AACCTTGG','ACACTGTGA'))
    print(output_lcs(lcs_backtrack(str1, str2), str1, str2))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
