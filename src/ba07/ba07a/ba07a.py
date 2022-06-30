"""
Rosalind: BA7A
Compute Distances Between Leaves

In this chapter, we define the length of a path in a tree as the sum of the lengths of its edges
(rather than the number of edges on the path).
As a result, the evolutionary distance between two present-day species
corresponding to leaves i and j in a tree T
is equal to the length of the unique path connecting i and j,
denoted d(i,j)(T).

Distance Between Leaves Problem
Compute the distances between leaves in a weighted tree.

Given: An integer n followed by the adjacency list of a weighted tree with n leaves.

Return: A space-separated nxn d(i,j),
where d(i,j) is the length of the path between leaves i and j.

Sample Dataset
4
0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6

Sample Output
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0

═════════════════════════════════════════════════

    [ Where am I? ]

Info.
  ┌─────────────────────────────────────────────────────────────────────────┐
  │ Distances Between Leaves Problem: (BA7A)                                │
  │ Compute the distances between leaves in a weighted tree.                │
  │   Input: A weighted tree with n leaves.                                 │
  │   Output: An nxn matrix (d(i,j)),                                       │
  │           where d(i,j) is the length of the path between leaves i and j.│
  └─────────────────────────────────────────────────────────────────────────┘

  * disstance matrix (fig 7.3)
    SPECIES   ALIGNMENT   DISTANCE MATRIX (Hamming distance)
                          Chimp   Human   Seal    Whale
    ---------------------------------------------------
    Chimp     ACGTAGGCCT  0       3       6       4
    Human     ATGTAAGACT  3       0       7       5
    Seal      TCGAGAGCAC  6       7       0       2
    Whale     TCGAAAGCAT  4       5       2       0

  * distance:
    a. the number of different characters between two strings
               v   v v   Hamming distance == 3
      Chimp - ACGTAGGCCT
      Human - ATGTAAGACT
    b. length of a path in a Tree, teh sum of the lengths of its edges (with weights)
       D(i,j) distance between leaves i and j
    c. D: distance matrix

  * fit:
    when d(i,j)(T) == D(i,j)
    T: weighted unrooted tree, D: distance matrix

  * addtive, non-addtive
    - additive: if there exists a tree that fits D, d(i,j)(T) == D(i,j)
    - non-additive: otherwise

  * Tree(D), simple tree:
    - tree that doesn't have a node with degree 2
    - Every simple Tree with n leaves has at most n-2 internal nodes.

  * searching for a Tree fitting a distance matrix
    1. three leaves

        3          D(1,3) = d(1,c) + d(3,c)
        |          D(2,3) = d(2,c) + d(3,c)
        c          D(1,2) = d(1,c) + d(2,c)
      /   \        internal node: c, leaves: 1,2,3
     1     2

     d(1,c) + d(2,c) = D(1,2)    --- (1)
     d(1,c) + d(3,c) = D(1,3)    --- (2)
     d(2,c) + d(3,c) = D(2,3)    --- (3)

     2d(1,c) + d(2,c) + d(3,c) = D(1,2) + D(1,3)      --- (1) + (2)
     2d(1,c) = D(1,2) + D(1,3) - (d(2,c) + d(3,c))
             = D(1,2) + D(1,3) - D(2,3)               --- by (3)
     d(1,c) = (D(1,2) + D(1,3) - D(2,3)) / 2          --- (4)

     2d(3,c) + d(1,c) + d(2,c) = D(1,3) + D(2,3)      --- (2) + (3)
     2d(3,c) = D(1,3) + D(2,3) - (d(1,c) + d(2,c))
             = D(1,3) + D(2,3) - D(1,2)               --- by (1)
     d(3,c) = (D(1,3) + D(2,3) - D(1,2)) / 2          --- (5)

     2d(2,c) + d(1,c) + d(3,c) = D(1,2) + D(2,3)      --- (1) + (3)
     2d(2,c) = D(1,2) + D(2,3) - (d(1,c) + d(3,c))
             = D(1,2) + D(2,3) - D(1,3)               --- by (2)
     d(2,c) = (D(1,2) + D(1,3) - D(1,3)) / 2          --- (6)

    2. four leaves
          i   j   k   l            i   j   k   l
       -----------------        ----------------
       i  0   13  21  22        i  0   3   4   3
       j  13  0   12  13        j  3   0   4   5
       k  21  12  0   13        k  4   4   0   2
       l  22  13  13  0         l  3   5   2   0
       additive                 non-additive (for this matrix)

                      D(i,j) =          d(i,a) + d(j,a) = 3
      i \     / k     D(i,k) = d(i,a) + d(a,b) + d(k,b) = 4
         a - b        D(i,l) = d(i,a) + d(a,b) + d(l,b) = 3
      j /     \ l     D(j,k) = d(j,a) + d(a,b) + d(k,b) = 4
                      D(j,l) = d(j,a) + d(a,b) + d(l,b) = 5
                      D(k,l) =          d(k,b) + d(l,b) = 2

                      D(i,j) = d(i,a) + d(a,b) + d(j,a) = 3
      i \     / j     D(i,k) =          d(i,a) + d(k,a) = 4
         a - b        D(i,l) = d(i,a) + d(a,b) + d(l,b) = 3
      k /     \ l     D(j,k) = d(j,b) + d(a,b) + d(k,a) = 4
                      D(j,l) =          d(j,b) + d(l,b) = 5
                      D(k,l) = d(k,a) + d(a,b) + d(l,b) = 2

                      D(i,j) = d(i,a) + d(a,b) + d(j,b) = 3
      i \     / j     D(i,k) = d(i,a) + d(a,b) + d(k,b) = 4
         a - b        D(i,l) =          d(i,a) + d(l,a) = 3
      l /     \ k     D(j,k) =          d(j,b) + d(k,b) = 4
                      D(j,l) = d(j,b) + d(a,b) + d(l,a) = 5
                      D(k,l) = d(k,b) + d(a,b) + d(l,a) = 2

           j          D(i,j) = d(i,c) + d(j,c) = 3
           |          D(i,k) = d(i,c) + d(k,c) = 4
       i - c - k      D(i,l) = d(i,c) + d(l,c) = 3
           |          D(j,k) = d(j,c) + d(k,c) = 4
           l          D(j,l) = d(j,c) + d(l,c) = 5
                      D(k,l) = d(k,c) + d(l,c) = 2

      D(i,j)+D(k,l), D(i,k)+D(j,l), D(i,l)+D(j,k) = 5, 9, 7 <-- does not satisfy four point condition

  * four point condition
    - d(i,j) + d(k,l) <= d(i,k) + d(j,l) = d(i,l) + d(j,k)
    - In terms of an nxn distance matrix,
      we say that a quartet of indices (i,j,l,k) satisfies the four point condition
      if two of the following sums are equal,
      and the third sum is less than or equal to the ohter two sums:
      (1) D(i,j) + D(k,l)    (2) D(i,k) + D(j,l)    (3) D(i,l) + D(j,k)

  * four point theorem
    A distance matrix is additive
    if and only if the four point condition holds for every quartet (i,j,k,l) of indices of the matrix.
    D(i,j) + D(k,l)    D(i,k) + D(j,l)    D(i,l) + D(j,k)

  ┌─────────────────────────────────────────────────────────────┐
  │ Distance-Based Phylogeny Problem:                           │
  │ Reconstruct an evolutionary tree fitting a distance matrix. │
  │   Input: A distance matrix.                                 │
  │   Output: A tree fitting this distance matrix.              │
  └─────────────────────────────────────────────────────────────┘
  TODO: requisites?

  * sample date

          1     2
          ↕     ↕
    0 <-> 4 <-> 5 <-> 3

  * What is the best representation of Weighted Adjacency List, Tree?

Plan 1.
  * graph to matrix, directly --> fill empty boxes
        0   1   2   3   4   5
      ┌───┬───┬───┬───┬───┬───┐      0->4:11    start   : 0,1,   2,3
    0 │ 0 │   │   │   │11 │   │      1->4:2     transfer: 4,     5
      ├───┼───┼───┼───┼───┼───┤      2->5:6     end     : 0,1,5, 2,3,4
    1 │   │ 0 │   │   │ 2 │   │      3->5:7
      ├───┼───┼───┼───┼───┼───┤      4->0:11    4->[0,1,5]
    2 │   │   │ 0 │   │   │ 6 │      4->1:2
      ├───┼───┼───┼───┼───┼───┤      4->5:4
    3 │   │   │   │ 0 │   │ 7 │      5->4:4     5->[2,3,4]
      ├───┼───┼───┼───┼───┼───┤      5->3:7
    4 │11 │ 2 │   │   │ 0 │ 4 │      5->2:6
      ├───┼───┼───┼───┼───┼───┤
    5 │   │   │ 6 │ 7 │ 4 │ 0 │      This is bidirectional directed graph.
      └───┴───┴───┴───┴───┴───┘      Two nodes point to each other.
                                     TODO: Is this a tree data structure?

    start -> transfer -> end
    0 ─┬─ 4 ─┬─ 0 -> (0,0) (1,0)
    1 ─┘     ├─ 1 -> (0,1) (1,1)
             └─ 5 -> (0,5) (1,5)
    2 ─┬─ 5 ─┬─ 2 -> (2,2) (3,2)
    3 ─┘     ├─ 3 -> (2,3) (3,3)
             └─ 4 -> (1,4) (3,4)

  * {(source,target):weight}
    {(0,4):11,(1,4):2,(2,5):6,(3,5):7,(4,0):11,(4,1):2,(4,5):4,(5,4):4,(5,3):7,(5,2):6}
  * {source:[(target,weight)]}
    {0:[(4,11)],1:[(4,2)],2:[(5,6)],3:[(5,7)],4:[(0,11),(1,2),(5,4)],5:[(4,4),(3,7),(2,6)]}

    d(0,0) = 0
    d(0,1) = d(0,4) + ┌ d(4,0)
                      ├ d(4,1)
                      └ d(4,5)
    d(0,2) = d(0,4) + ┌ d(4,0)
                      ├ d(4,1) + d(1,4)
                      └ d(4,5) + ┌ d(5,2)
                                 ├ d(5,3)
                                 └ d(5,4)

    d(i,j) = ┌ if i == j                   0                             e.g. d(i,i)
             ├ if (i,j) in dictionary      d(i,j)                        e.g. d(0,4)
             ├ if len(d(i,_)) == 1 OR len(d(_,j)) == 1   (just one transfer station)
             │  ┌ if len((i,_)) == 1       d(i,i1) + d(i1,j)             e.g. d(2,3) = d(2,5) + d(5,3)
             │  └ if len((i,_)) != 1       d(i,i1) + d(i1,j1) + d(j1,j)  e.g. d(3,1) = d(3,5) + d(5,4) + d(4,1)
             └ if len(d(i,_)) != 1 AND len(d(_,j)) != 1  (multiple transfer stations)
                 d(i,j) = d(i,_) + d(_,_) + ... + d(_,_) + d(_,j)      complex and difficult

  * TODO: Can I solve the problems of the remaining chapters with a simple dictionary?
          Should I create some abstract data structures for following quizes?

═════════════════════════════════════════════════

References:
- Split string with multiple delimiters in Python [duplicate]
  https://stackoverflow.com/questions/4998629/split-string-with-multiple-delimiters-in-python
- Creating a 2d matrix in python
  https://stackoverflow.com/questions/4230000/creating-a-2d-matrix-in-python
- Find all paths in DAG from starting node **
  https://stackoverflow.com/questions/54441582/find-all-paths-in-dag-from-starting-node
- Toward an Algorithm for Distance-Based Phylogeny Construction
  https://www.youtube.com/watch?v=TSIkxBoaykA
"""
#!/usr/bin/env python
import time
import re


#################################################
# dictionary
#################################################


def allpaths_old(graph, start, sink, path=[]):
    """
    ({a:[a]},a,a,[a]) -> [[a]]
    returns all path from 'start' node to 'sink' node (BA5D)
    https://stackoverflow.com/questions/54441582/find-all-paths-in-dag-from-starting-node
    This is very useful function.
    >>> allpaths_old({'0':['1','2'],'1':['4'],'2':['3'],'3':['4'],'5':['1']},0,4)
        [['0','1','4'],['0','2','3','4']]
    >>> allpaths_old({0:[4],1:[4],2:[5],3:[5],4:[0,1,5],5:[4,3,2]})
        []  # <-- It doesn't work in bidirectional directed graph.
    """
    path = path + [start]
    if start not in graph:
        return [path]
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = allpaths(graph, node, sink, path)
            for newpath in newpaths:
                if newpath[-1] == sink:
                    paths.append(newpath)
    return paths


def allpaths(graph, start, sink, path=[]):
    """
    ({a:[a]},a,a,[a]) -> [[a]]
    returns all path from 'start' node to 'sink' node (BA5D)
    modified version to solve BA7A
    >>> allpaths({'0':['1','2'],'1':['4'],'2':['3'],'3':['4'],'5':['1']},'0','4')
        [['0','1','4'],['0','2','3','4']]
    >>> allpaths({0:[4],1:[4],2:[5],3:[5],4:[0,1,5],5:[4,3,2]},0,3)
        [[0,4,5,3]]
    """
    path = path + [start]
    if start not in graph:
        return [path]
    paths = []
    for node in graph[start]:
        if path[-1] == sink and path not in paths:
            paths.append(path)
        if node not in path:
            newpaths = allpaths(graph, node, sink, path)
            for newpath in newpaths:
                if newpath[-1] == sink:
                    paths.append(newpath)
    return paths


def distance_matrix(n, graph, graph_wt):
    """
    (int,{a:[a]},{(a,a):int})
    returns a matrix of nxn d(i,j) values
    """
    matrix = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            path = allpaths(graph,i,j)[0]  # [[0, 4, 5, 3]]
            distance = 0
            for k in range(len(path)-1):
                distance += graph_wt[(path[k],path[k+1])]
            matrix[i][j] = distance
    return matrix


#################################################
# TODO: OOP version
#################################################


def main():
    f = open('/home/wsl/rosalind/data/ba07a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    graph_wt = {}
    graph = {}
    n = int(lines[0])
    for line in lines[1:]:
        source, target, weight = [int(elem) for elem in re.split('->|:',line)]
        graph_wt[(source,target)] = weight
        if source not in graph:
            graph[source] = [target]
        else:
            graph[source].append(target)

    start_time = time.time()
    matrix = distance_matrix(n, graph, graph_wt)
    for row in matrix:
        print(' '.join(map(str, row)))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()