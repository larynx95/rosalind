"""
Rosalind: BA5J
Align Two Strings Using Affine Gap Penalties

A gap is a contiguous sequence of spaces in a row of an alignment.
One way to score gaps more appropriately is to define an affine penalty for a gap of length k as σ + ε · (k - 1),
where σ is the gap opening penalty, assessed to the first symbol in the gap,
and ε is the gap extension penalty, assessed to each additional symbol in the gap.
We typically select ε to be smaller than σ so that the affine penalty
for a gap of length k is smaller than the penalty for k independent single-nucleotide indels (σ · k).

Alignment with Affine Gap Penalties Problem
Construct a highest-scoring global alignment of two strings (with affine gap penalties).

Given: Two amino acid strings v and w (each of length at most 100).

Return:
The maximum alignment score between v and w,
followed by an alignment of v and w achieving this maximum score.
Use the "BLOSUM62" scoring matrix,
a gap opening penalty(σ) of 11,
and a gap extension penalty(ε) of 1.


Sample Dataset
PRTEINS
PRTWPSEIN

Sample Output
8
PRT---EINS
PRTWPSEIN-

═════════════════════════════════════════════════

Info.
  * What is the "gap", "affine penalty"?
    - affine penalty (gap length k): sigma + epsilon * (k - 1)
    - example
      GATCCAG                    GATCCAG
      GA-C-AG                    GA--CAG
        ^ ^                        ^^
        │ └ gap opening            │└ gap extension (epsilon, ε)
        └── gap opening            └─ gap opening (sigma, σ)

      2*(σ + ε*(1-1))            σ + ε*(2-1)
      = 2σ = 10                  = 5 + 1 = 6      if σ=5, ε=1

  * sample dataset
    8
    PRT---EINS
    PRTWPSEIN-

Plan 1.

═════════════════════════════════════════════════

References:
- Week 2.2.2 Alignment with Affine Gap Penalty and Calculation of Time Complexity of The Needleman Wun
  https://www.youtube.com/watch?v=DQQ_q2dn2ds
- Gap Penalties CMSC 423
  https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
- Bioinformatics 1: lecture 4
  https://www.cs.purdue.edu/homes/ayg/TALKS/BLORE10/lecture4.pdf
- Advanced Sequence Alignment
  http://www.csbio.unc.edu/mcmillan/Comp555S16/Lecture14.html
"""

import time


dic_blosum62 = {'A': [ 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2],
                'C': [ 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
                'D': [-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3],
                'E': [-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2],
                'F': [-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3],
                'G': [ 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3],
                'H': [-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2],
                'I': [-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1],
                'K': [-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2],
                'L': [-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1],
                'M': [-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1],
                'N': [-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2],
                'P': [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3],
                'Q': [-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1],
                'R': [-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2],
                'S': [ 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2],
                'T': [ 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2],
                'V': [ 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1],
                'W': [-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2],
                'Y': [-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7]}


def main():
    f = open('/home/wsl/rosalind/data/ba05j.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
