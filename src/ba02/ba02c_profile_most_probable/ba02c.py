"""
Rosalind: BA2C
Find a Profile-most Probable k-mer in a String

Given a profile matrix Profile,
we can evaluate the probability of every k-mer in a string Text
and find a Profile-most probable k-mer in Text,
i.e., a k-mer that was most likely to have been generated
by Profile among all k-mers in Text.
For example, ACGGGGATTACC is the Profile-most probable 12-mer in GGTACGGGGATTACCT.
Indeed, every other 12-mer in this string has probability 0.

In general, if there are multiple Profile-most probable k-mers in Text,
then we select the first such k-mer occurring in Text.

Given: A string Text, an integer k, and a 4 * k matrix Profile.

Return: A Profile-most probable k-mer in Text.
(If multiple answers exist, you may return any one.)

Sample Dataset
ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
5
0.2 0.2 0.3 0.2 0.3
0.4 0.3 0.1 0.5 0.1
0.3 0.3 0.5 0.2 0.4
0.1 0.2 0.1 0.1 0.2

Sample Output
CCGAG

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference? (very important)
      -> PREV: Median string problem (BA2B) - not fast enough... find another algorithm
      ↓
    * GreedyMotifSearch
      -> HERE: Profile-most probable k-mer problem (BA2C)
      -> NEXT: GreedyMotifSearch (BA2D)

Info:
    Profile-most Probable k-mer Problem
    Find a Profile-most probable k-mer in a string.
    (profile-most probable k-mer is not consensus.)

                          T  C  G  G  G  G  g  T  T  T  t  t
                          c  C  G  G  t  G  A  c  T  T  a  C
                          a  C  G  G  G  G  A  T  T  T  t  C
                          T  t  G  G  G  G  A  c  T  T  t  t
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C
                          T  t  G  G  G  G  A  c  T  T  C  C
                          T  C  G  G  G  G  A  T  T  c  a  t
                          T  C  G  G  G  G  A  T  T  c  C  t
                          T  a  G  G  G  G  A  a  c  T  a  C
                          T  C  G  G  G  t  A  T  a  a  C  C

    SCORE(Motifs)         3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30

    COUNT(Motifs)     A:  2  2  0  0  0  0  9  1  1  1  3  0
                      C:  1  6  0  0  0  0  0  4  1  2  4  6
                      G:  0  0 10 10  9  9  1  0  0  0  0  0
                      T:  7  2  0  0  1  1  0  5  8  7  3  4

    PROFILE(Motifs)   A: .2 .2  0  0  0  0 .9 .1 .1 .1 .3  0
                      C: .1 .6  0  0  0  0  0 .4 .1 .2 .4 .6
                      G:  0  0  1  1 .9 .9 .1  0  0  0  0  0
                      T: .7 .2  0  0 .1 .1  0 .5 .8 .7 .3 .4
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

                      A: .2 .2 .0 .0 .0 .0 .9 .1 .1 .1 .3 .0
                      C: .1 .6 .0 .0 .0 .0 .0 .4 .1 .2 .4 .6
                      G: .0 .0  1  1 .9 .9 .1 .0 .0 .0 .0 .0
                      T: .7 .2 .0 .0 .1 .1 .0 .5 .8 .7 .3 .4
    Pr(ACGGGGATTACC|Profile) = .2*.6*1*1*.9*.9*.9*.5*.8*.1*.4*.6 = 0.000839808
    (profile most probable kmer is not the same as concensus!!)

Plan 1.
- What is the "Greedy algorithm"?
  ; A greedy algorithm is an approach for solving a problem
    by selecting the best option available at the moment.
    It doesn't worry whether the current best result will bring the overall optimal result.
- Profile-most probable k-mer finding is a kind of greedy algorithm.
"""

import time


def all_kmers(dna, k):
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


def pr_score(profile, kmer):
    """
    ([[float]],str) -> float
    profile most probable kmer
    >>> profile = [[0.2,0.2,0.3,0.2,0.3],[0.4,0.3,0.1,0.5,0.1],[0.3,0.3,0.5,0.2,0.4],[0.1,0.2,0.1,0.1,0.2]]
    >>> pr_score(profile,"ACGTC")
        0.00030000000000000003
    """
    val = 1
    for i in range(len(kmer)):
        nuc = kmer[i]
        if nuc =='A':
            val *= profile[0][i]
        elif nuc == 'C':
            val *= profile[1][i]
        elif nuc == 'G':
            val *= profile[2][i]
        elif nuc =='T':
            val *= profile[3][i]
    return val


def prof_most_probable_kmer(dna, k, profile):
    """
    (str,int,[[float]]) -> {str}
    get profile most probable k-mer from a DNA string and given profile
    ignoring the fact that there can be more than one k-mer
    TODO: Improve this function
    TODO: Implement "Motifs(profile, DNAs)" function later
    >>> profile = [[0.2,0.2,0.3,0.2,0.3],[0.4,0.3,0.1,0.5,0.1],[0.3,0.3,0.5,0.2,0.4],[0.1,0.2,0.1,0.1,0.2]]
    >>> dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
    >>> prof_most_probable_kmer(dna,5,profile)
        {'CCGAG'}
    """
    dict = {}
    largest = 0
    for kmer in all_kmers(dna, k):
        val = pr_score(profile, kmer)
        dict[kmer] = val
        if val > largest:
            largest = val
    return {k for k, v in dict.items() if v == largest}


def main():
    # this part is just for reading a file and getting data
    f = open('/home/wsl/rosalind/data/ba02c.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    dna = lines[0]
    k = int(lines[1])
    profile = []
    for line in lines[2:]:
        profile.append([float(num) for num in line.strip().split()])
    f.close()

    # display result
    start_time = time.time()
    print(prof_most_probable_kmer(dna, k, profile))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()