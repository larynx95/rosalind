"""
Rosalind: BA5B
Find the Length of a Longest Path in a Manhattan-like Grid

Length of a Longest Path in the Manhattan Tourist Problem
Find the length of a longest path in a rectangular city.

Given:
Integers n and m,
followed by an n * (m+1) matrix Down
and an (n+1) * m matrix Right.
The two matrices are separated by the "-" symbol.

Return:
The length of a longest path from source (0,0) to sink (n,m)
in the n * m rectangular grid
whose edges are defined by the matrices Down and Right.

Sample Dataset
4 4
1 0 2 4 3
4 6 5 2 1
4 4 5 2 1
5 6 8 5 3
-
3 2 4 0
3 2 4 2
0 7 3 3
3 3 0 2
1 3 2 2

Sample Output
34

═════════════════════════════════════════════════

    [ Where am I? ]

    * PREV: Dynamic programming (DP) in "coin change problem" (BA5A)
      ↓
    * HERE: Applying dynamic programming to "Manhattan tourist problem" - Longest Path (BA5B)
      ↓
    * "Manhattan Tourist Problem" to "Longest Common Subsequence (LCS) Problem"
      -> LCS problem (BA5C) & Edit distance of two strings (BA5G)
      -> Generalized version of algorithm for any Directed Acyclic Graph (DAG) (BA5D)
         -> Topological ordering (BA5N)
      ↓
    * Alignment of Two Strings Problems
      -> Global alignment (BA5E)
         -> Affine alingment (BA5J)
      -> Local alignment (BA5F)
         -> Fitting alignment (BA5H)
         -> Overlap alignment (BA5I)



Plan 1.
  * What is the n * m rectangular grid?
      ┌─┬─┬─┐ Is this 4*4 grid or 3*3 grid? 3*3 grid
      ├─┼─┼─┤ n*m grid has n+1 vertical lines, m+1 horizontal lines.
      ├─┼─┼─┤
      └─┴─┴─┘
  * sample dataset (4*4 grid):
    (4+1 vertical lines, 4+1 horizonatal lines)  ← ↑ → ↓ ↖ ↗ ↘ ↙
      ┌──3──┬──2──┬──4──┬──0──┐    0─────3─────5─────9─────9    0─────3─────5─────9─────9
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      1     0     2     4     3    ^     │     │     │     │    v     │     │     │     │
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      ├──3──┼──2──┼──4──┼──2──┤    1 <-  4─────7────13────15    1 ->  4─────7────13────15
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      4     6     5     2     1    │     ^     │     │     │    │     v     │     │     │
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      ├──0──┼──7──┼──3──┼──3──┤    5────10 <- 17────20────23    5────10 -> 17────20────23
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      4     4     5     2     1    │     │     ^     │     │    │     │     v     |     │
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      ├──3──┼──3──┼──0──┼──2──┤    9────14────22────22────24    9────14────22────22────24
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      5     6     8     5     3    │     │     ^     │     │    │     │     v     │     │
      │     │     │     │     │    │     │     │     │     │    │     │     │     │     │
      └──1──┴──3──┴──2──┴──2──┘    14───20────30 <- 32 <- 34    14───20────30 -> 32 -> 34
      matirx                       barcktracking                longest path

      dynamic table
      [[0,   3,  5,  9,  9],
       [1,   4,  7, 13, 15],
       [5,  10, 17, 20, 23],
       [9,  14, 22, 22, 24],
       [14, 20, 30, 32, 34]]

  * expressed as (i,j) nodes: ← ↑ → ↓ ↖ ↗ ↘ ↙
      00 - 10 - 20 - 30 - 40
      |    |    |    |    |
      01 - 11 - 21 - 31 - 41
      |    |    |    |    |
      02 - 12 - 22 - 32 - 42
      |    |    |    |    |
      03 - 13 - 23 - 33 - 43
      |    |    |    |    |
      04 - 14 - 24 - 34 - 44
      |    |    |    |    |
      05 - 15 - 25 - 35 - 45

  * algorithms
  - recursive 'SOUTHOREAST' algorithm suffers from a huge number of recursive calls
  ╔═════════════════════════════════════════════════════════════════════════╗
  ║ SOUTHOREAST(i, j)                                                       ║
  ║     if i = 0 and j = 0                                                  ║
  ║         return 0                                                        ║
  ║     x <- -infinite, y <- -infinite                                      ║
  ║     if i > 0                                                            ║
  ║         x <- SOUTHOREAST(i-1,j) + weight of vertical edge into (i,j)    ║
  ║     if j > 0                                                            ║
  ║         y <- SOUTHOREAST(i,j-1) + weight of horizontal edge into (i,j)  ║
  ║     return max{x, y}                                                    ║
  ╚═════════════════════════════════════════════════════════════════════════╝

  ╔═════════════════════════════════════════════════════════════════════════╗
  ║ MANHATTANTOURIST(n, m, Down, Right)                                     ║
  ║     s(0,0) <- 0                                                         ║
  ║     for i <- 1 to n                                                     ║
  ║         s(i,0) <- s(i-1,0) + down(i,0)                                  ║
  ║     for j <- 1 to m                                                     ║
  ║         s(0,j) <- s(0,j-1) + right(0,j)                                 ║
  ║     for i <- 1 to n                                                     ║
  ║         for j <- 1 to m                                                 ║
  ║             s(i,j) <- max{s(i-1,j) + down(i,j), s(i,j-1) + right(i,j)}  ║
  ║     return s(n,m)                                                       ║
  ╚═════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
"""

import time


def south_or_east(row, col, down, right):
    """
    (int,int,[[int]],[[int]])
    returns the length of a longest path from source (0,0) to sink (row, col)
    This recursive algorithm suffers from a huge number of recursive calls
    >>> south_or_east(2,2,[[1,0,2,4,3],[4,6,5,2,1],[4,4,5,2,1],[5,6,8,5,3]],
                          [[3,2,4,0],[3,2,4,2],[0,7,3,3],[3,3,0,2],[1,3,2,2]])
        17
    """
    if row == 0 and col == 0:
        return 0
    x, y = -float('inf'), -float('inf')
    if row > 0:
        x = south_or_east(row-1, col, down, right) + down[row-1][col]
    if col > 0:
        y = south_or_east(row, col-1, down, right) + right[row][col-1]
    return max(x,y)


def manhattan_tourist(row, col, down, right):
    """
    (int,int,[[int]],[[int]]) -> int
    returns longest path
    algorithm in textbook, using dynamic programming, rectangular grid
    >>> manhattan_tourist(4,4,[[1,0,2,4,3],[4,6,5,2,1],[4,4,5,2,1],[5,6,8,5,3]],
                              [[3,2,4,0],[3,2,4,2],[0,7,3,3],[3,3,0,2],[1,3,2,2]])
        34
    """
    nodes = [[0] * (col+1) for i in range(row+1)]   # dynamic table
    for i in range(1,row+1):                        # add all values of vertical edges
        nodes[i][0] = nodes[i-1][0] + down[i-1][0]
    for j in range(1,col+1):                        # add all values of horizontal edges
        nodes[0][j] = nodes[0][j-1] + right[0][j-1]
    for i in range(1,row+1):                        # add all values of other edges
        for j in range(1,col+1):
            nodes[i][j] = max(nodes[i-1][j] + down[i-1][j], nodes[i][j-1] + right[i][j-1])
    return nodes[row][col]                          # return a value of a specific node




def main():
    # read data
    try:
        with open('/home/wsl/rosalind/data/ba05b.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            n, m = [int(line.strip()) for line in lines[0].split(' ')]
            vlines = []
            nline = 1
            for line in lines[1:]:
                nline += 1
                if line == '-':
                    break
                vlines.append([int(elem) for elem in line.split(' ')])
            hlines = []
            for line in lines[nline:]:
                if line == ' ' or line == '\n':
                    break
                hlines.append([int(elem) for elem in line.split(' ')])
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    # execute code
    start_time = time.time()
    print(manhattan_tourist(4,4,[[1,0,2,4,3],[4,6,5,2,1],[4,4,5,2,1],[5,6,8,5,3]],\
                                [[3,2,4,0],[3,2,4,2],[0,7,3,3],[3,3,0,2],[1,3,2,2]]))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
