"""
Rosalind: BA5L
Align Two Strings Using Linear Space

The pseudocode below for LinearSpaceAlignment describes
how to recursively find a longest path in the alignment graph constructed
for a substring vtop+1 ... vbottom of v and a substring wleft+1 ... wright of w.
LinearSpaceAlignment calls the function MiddleNode(top, bottom, left, right),
which returns the coordinate i of the middle node (i, j) defined by the sequences vtop+1 ... vbottom and wleft+1 ... wright.
LinearSpaceAlignment also calls MiddleEdge(top, bottom, left, right), which returns → , ↓, or ↘
depending on whether the middle edge is horizontal, vertical, or diagonal.
The linear-space alignment of strings v and w is constructed by calling LinearSpaceAlignment(0, n, 0, m).
The case left = right describes the alignment of an empty string against the string vtop+1 ... vbottom,
which is trivially computed as the score of a gap formed by bottom - top vertical edges.

╔════════════════════════════════════════════════════════════════════╗
║ LinearSpaceAlignment(top, bottom, left, right)                     ║
║     if left = right                                                ║
║         return alignment formed by bottom - top vertical edges     ║
║     if top = bottom                                                ║
║         return alignment formed by right - left horizontal edges   ║
║     middle <- └ (left + right)/2 ┘                                 ║
║     midNode <- MiddleNode(top, bottom, left, right)                ║
║     midEdge <- MiddleEdge(top, bottom, left, right)                ║
║     LinearSpaceAlignment(top, midNode, left, middle)               ║
║     output midEdge                                                 ║
║     if midEdge = "→" or midEdge = "↘"                             ║
║         middle <- middle + 1                                       ║
║     if midEdge = "↓" or midEdge ="↘"                              ║
║         midNode <- midNode + 1                                     ║
║     LinearSpaceAlignment(midNode, bottom, middle, right)           ║
║                                                                    ║
╚════════════════════════════════════════════════════════════════════╝

Global Alignment in Linear Space Problem
Find the highest-scoring alignment between two strings using a scoring matrix in linear space.

Given: Two long amino acid strings (of length approximately 10,000).

Return: The maximum alignment score of these strings, followed by an alignment achieving this maximum score.
Use the BLOSUM62 scoring matrix and indel penalty σ = 5.

Sample Dataset
PLEASANTLY
MEANLY

Sample Output
8
PLEASANTLY
-MEA--N-LY

═════════════════════════════════════════════════

Plan 1.

═════════════════════════════════════════════════

References:
"""

import time


def main():
    f = open('/home/wsl/rosalind/data/ba05l.txt', 'r')
    lines = [line.strip() for line in f.readlines()]

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
