"""
Rosalind: BA7E
Implement the Neighbor Joining Algorithm

The pseudocode below summarizes the neighbor-joining algorithm.

NeighborJoining(D,n)
    if n = 2
        T <- tree consisting of a single edge of length D1,2
        return T
    D' <- neighbor-joining matrix constructed from the distance matrix D
    find elements i and j such that D'i,j is a minimum non-diagonal element of D'
    Δ <- (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
    limbLengthi <- (1/2)(Di,j + Δ)
    limbLengthj <- (1/2)(Di,j - Δ)
    add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
    remove rows i and j from D
    remove columns i and j from D
    T <- NeighborJoining(D, n - 1)
    add two new limbs (connecting node m with leaves i and j) to the tree T
    assign length limbLengthi to Limb(i)
    assign length limbLengthj to Limb(j)
    return T

Neighbor Joining Problem
Construct the tree resulting from applying the neighbor-joining algorithm to a distance matrix.

Given: An integer n, followed by a space-separated n x n distance matrix.

Return: An adjacency list for the tree resulting from applying the neighbor-joining algorithm.
Edge-weights should be accurate to two decimal places
(they are provided to three decimal places in the sample output below).

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
The n leaves must be labeled 0, 1, ..., n-1 in order of their appearance in the distance matrix.
Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.

Sample Dataset
4
0   23  27  20
23  0   30  28
27  30  0   30
20  28  30  0

Sample Output
0->4:8.000
1->5:13.500
2->5:16.500
3->4:12.000
4->5:2.000
4->0:8.000
4->3:12.000
5->1:13.500
5->2:16.500
5->4:2.000

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
    f = open('/home/wsl/rosalind/data/ba07e.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()