"""
Rosalind: BA7G
Adapt SmallParsimony to Unrooted Trees

When the position of the root in an evolutionary tree is unknown,
we can simply assign the root to any edge that we like,
apply SmallParsimony from “Implement SmallParsimony” to the resulting rooted tree, and then remove the root.
It can be shown that this method provides a solution to the following problem.

Small Parsimony in an Unrooted Tree Problem
Find the most parsimonious labeling of the internal nodes in an unrooted tree.

Given: An unrooted binary tree with each leaf labeled by a string of length m.

Return: A labeling of all other nodes of the tree by strings of length m
that minimizes the parsimony score of the tree.

Note on formatting: Your internal node labelings may differ from the sample provided.

Sample Dataset
4
TCGGCCAA->4
4->TCGGCCAA
CCTGGCTG->4
4->CCTGGCTG
CACAGGAT->5
5->CACAGGAT
TGAGTACC->5
5->TGAGTACC
4->5
5->4

Sample Output
17
TCGGCCAA->CCAGGCAC:4
CCTGGCTG->CCAGGCAC:3
TGAGTACC->CAAGGAAC:4
CCAGGCAC->CCTGGCTG:3
CCAGGCAC->CAAGGAAC:2
CCAGGCAC->TCGGCCAA:4
CACAGGAT->CAAGGAAC:4
CAAGGAAC->CACAGGAT:4
CAAGGAAC->TGAGTACC:4
CAAGGAAC->CCAGGCAC:2

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
    f = open('/home/wsl/rosalind/data/ba07g.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()