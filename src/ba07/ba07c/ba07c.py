"""
Rosalind: BA7C
Implement AdditivePhylogeny

The following recursive algorithm, called AdditivePhylogeny,
finds the simple tree fitting an n x n additive distance matrix D.
We assume that you have already implemented a program Limb(D, j)
that computes LimbLength(j) for a leaf j based on the distance matrix D.
Rather than selecting an arbitrary leaf j from Tree(D) for trimming,
AdditivePhylogeny selects leaf n (corresponding to the last row and column of D).

╔════════════════════════════════════════════════════════════════════════════════════════════╗
║ AdditivePhylogeny(D, n)                                                                    ║
║     if n = 2                                                                               ║
║         return the tree consisting of a single edge of length D1,2                         ║
║     limbLength <- Limb(D, n)                                                               ║
║     for j <- 1 to n - 1                                                                    ║
║         Dj,n <- Dj,n - limbLength                                                          ║
║         Dn,j <- Dj,n                                                                       ║
║     (i,n,k) <- three leaves such that Di,k = Di,n + Dn,k                                   ║
║     x <- Di,n                                                                              ║
║     remove row n and column n from D                                                       ║
║     T <- AdditivePhylogeny(D, n - 1)                                                       ║
║     v <- the (potentially new) node in T at distance x from i on the path between i and k  ║
║     add leaf n back to T by creating a limb (v, n) of length limbLength                    ║
║     return T                                                                               ║
╚════════════════════════════════════════════════════════════════════════════════════════════╝

Additive Phylogeny Problem
Construct the simple tree fitting an additive matrix.

Given: n and a tab-delimited n x n additive matrix.

Return: A weighted adjacency list for the simple tree fitting this matrix.

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
The n leaves must be labeled 0, 1, ..., n-1 in order of their appearance in the distance matrix.
Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.

Sample Dataset
4
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0

Sample Output
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
    f = open('/home/wsl/rosalind/data/ba07c.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()