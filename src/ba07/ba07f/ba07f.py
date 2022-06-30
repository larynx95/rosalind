"""
Rosalind: BA7F
Implement SmallParsimony

The pseudocode for SmallParsimony is shown below.
It returns the parsimony score for a binary rooted tree T
whose leaves are labeled by symbols stored in an array Character (i.e., Character(v) is the label of leaf v).
At each iteration, it selects a node v and computes sk(v) for each symbol k in the alphabet.
For each node v, SmallParsimony maintains a value Tag(v),
which indicates whether the node has been processed
(i.e., Tag(v) = 1 if the array sk(v) has been computed and Tag(v) = 0 otherwise).
We call an internal node of T ripe if its tag is 0 but its children’s tags are both 1.
SmallParsimony works upward from the leaves, finding a ripe node v at which to compute sk(v) at each step.

SmallParsimony(T, Character)
    for each node v in tree T
        Tag(v) <- 0
        if v is a leaf
            Tag(v) <- 1
            for each symbol k in the alphabet
                if Character(v) = k
                    sk(v) <- 0
                else
                    sk(v) <- ∞
    while there exist ripe nodes in T
        v <- a ripe node in T
        Tag(v) <- 1
        for each symbol k in the alphabet
            sk(v) <- minimum over all symbols i {si(Daughter(v))+δi,k} + minimum over all symbols j {sj(Son(v))+δj,k}
   return minimum over all symbols k {sk(v)}

Small Parsimony Problem
Find the most parsimonious labeling of the internal nodes of a rooted tree.

Given: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.

Return: The minimum parsimony score of this tree,
followed by the adjacency list of the tree corresponding to labeling internal nodes by DNA strings
in order to minimize the parsimony score of the tree.

Note: Remember to run SmallParsimony on each individual index of the strings at the leaves of the tree.

Sample Dataset
4
4->CAAATCCC
4->ATTGCGAC
5->CTGCGCTG
5->ATGGACGA
6->4
6->5

Sample Output
16
ATTGCGAC->ATAGCCAC:2
ATAGACAA->ATAGCCAC:2
ATAGACAA->ATGGACTA:2
ATGGACGA->ATGGACTA:1
CTGCGCTG->ATGGACTA:4
ATGGACTA->CTGCGCTG:4
ATGGACTA->ATGGACGA:1
ATGGACTA->ATAGACAA:2
ATAGCCAC->CAAATCCC:5
ATAGCCAC->ATTGCGAC:2
ATAGCCAC->ATAGACAA:2
CAAATCCC->ATAGCCAC:5

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
    f = open('/home/wsl/rosalind/data/ba07f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()