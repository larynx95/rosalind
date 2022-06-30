"""
Rosalind: BA7D
Implement UPGMA

The pseudocode for UPGMA is shown below.

UPGMA(D, n)
    Clusters ← n single-element clusters labeled 1, ... , n
    construct a graph T with n isolated nodes labeled by single elements 1, ... , n
    for every node v in T
        Age(v) ← 0
    while there is more than one cluster
        find the two closest clusters Ci and Cj
        merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        add a new node labeled by cluster Cnew to T
        connect node Cnew to Ci and Cj by directed edges
        remove the rows and columns of D corresponding to Ci and Cj
        remove Ci and Cj from Clusters
        add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        add Cnew to Clusters
    root ← the node in T corresponding to the remaining cluster
    for each edge (v, w) in T
        length of (v, w) ← Age(v) - Age(w)
    return T

UPGMA Problem
Construct the ultrametric tree resulting from UPGMA.

Given: An integer n followed by a space-delimited n x n distance matrix.

Return: An adjacency list for the ultrametric tree output by UPGMA.
Weights should be accurate to three decimal places.

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
The n leaves must be labeled 0, 1, ..., n-1 in order of their appearance in the distance matrix.
Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.

Sample Dataset
4
0   20  17  11
20  0   20  13
17  20  0   10
11  13  10  0

Sample Output
0->5:7.000
1->6:8.833
2->4:5.000
3->4:5.000
4->2:5.000
4->3:5.000
4->5:2.000
5->0:7.000
5->4:2.000
5->6:1.833
6->5:1.833
6->1:8.833

═════════════════════════════════════════════════

Info.

Plan 1.

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time


def main():
    f = open('/home/wsl/rosalind/data/ba07d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()