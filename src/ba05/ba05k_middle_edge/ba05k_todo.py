"""
Rosalind: BA5K
Find a Middle Edge in an Alignment Graph in Linear Space

Middle Edge in Linear Space Problem
Find a middle edge in the alignment graph in linear space.

Given:
Two amino acid strings.

Return:
A middle edge in the alignment graph of these strings,
where the optimal path is defined by the BLOSUM62 scoring matrix and a linear indel penalty equal to 5.
Return the middle edge in the form “(i, j) (k, l)”, where (i, j) connects to (k, l).

Sample Dataset
PLEASANTLY
MEASNLY

Sample Output
(4, 3) (5, 4)

═════════════════════════════════════════════════

Plan 1.

═════════════════════════════════════════════════

References:
"""

import time


def main():
    f = open('/home/wsl/rosalind/data/ba05k.txt', 'r')
    lines = [line.strip() for line in f.readlines()]

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
