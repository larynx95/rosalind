"""
The BEST theorem

1. Matrix representation of Eulerian Graph

  De Bruijn Graph (G):    {a:[b], b:[c,d], c:[a,b], d:[c]}

  (1) matrix A(G):
                   target
                   a  b  c  d
                a  0  1  0  0
      source    b  0  0  1  1
                c  1  1  0  0
                d  0  0  1  0
      sum(col)     1  2  2  1   <-- Indegree(node)

  (2) matrix -A(G):
                   target
                   a  b  c  d
                a  0 -1  0  0
      source    b  0  0 -1 -1
                c -1 -1  0  0
                d  0  0 -1  0

  (3) matrix A*(G):
                   target
                   a  b  c  d    replacing i-th diagonal entry of -A(G)
                a  1 -1  0  0    by Indegree(i) for each node i in G
      source    b  0  3 -1 -1    c(G) == 2
                c -1 -1  2  0
                d  0  0 -1  1

  (4) cofactor table

    ┌ 2  2  2  2 ┐  For a given Eulerian graph G, all i-cofactors of A*(G) have the same value.
    │ 2  2  2  2 │
    │ 2  2  2  2 │
    └ 2  2  2  2 ┘

  (5) the number of Eulerian cycles

    n = c(G) * ∏ <all nodes v in graph G> (InDegree(v) - 1)!

                   target
                   a  b  c  d
                a  0  1  0  0
      source    b  0  0  1  1
                c  1  1  0  0
                d  0  0  1  0
      sum(col)     1  2  2  1   <-- Indegree(node)

    n = 2 * sum ┌ v == a, (1-1)! = 2  So there are two Eulerian cycles.
                ├ v == b, (2-1)!
                ├ v == c, (2-1)!
                └ v == d, (1-1)!

2. i-cofactor

  i-cofactor:
    ┌ 00[01]02 03 ┐  ->  ┌ -  *  -  -  ┐  ->  ┌ 10 12 13 ┐ -> calculate determinant
    │ 10 11 12 13 │      │ 10 -  12 13 │      │ 20 22 23 │
    │ 20 21 22 23 │      │ 20 -  22 23 │      └ 30 32 33 ┘
    └ 30 31 32 33 ┘      └ 30 -  32 33 ┘

  determinant:
    ┌ a b c ┐ --> |A| = a(ei - fh) - b(di - fg) + c(dh - eg)
    │ d e f │
    └ g h i ┘

═════════════════════════════════════════════════

References:
- Cofactor: Definition & Formula
  https://study.com/academy/lesson/cofactor-definition-formula.html
"""

def graph_to_matrix(graph):
    """
    returns a matrix from De Gruijn Graph
    TODO:
    """
    pass


def determinant(matrix):
    """
    returns determinant of a matrix
    TODO:
    """
    pass


def cofactor(matrix, row, col):
    """
    returns a coractor table
    TODO:
    """
    pass


def num_eulerian_cycles(graph):
    """
    returns the number of Eulerian cycles
    TODO:
    """
    pass
