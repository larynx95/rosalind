"""
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

    MOTIFENUMERATION(Dna, k, d)
        Patterns <- an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern' differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns

Implanted Motif Problem
Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.

Given: Integers k and d, followed by a collection of strings Dna.

Return: All (k, d)-motifs in Dna.

Sample Dataset
3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output
ATA ATT GTT TTT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    * find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> HERE: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference?
      -> NEXT: "Median string problem" (BA2B)

Info:
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
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

Thoughts:
- What exactly is the meaning of this exercise? What are the problems?
  - finding motifs from multiple DNA string (e.g. 15-mer with mutations from each of ten DNA string )
  - Can this exercise be solved by functions in BA1B or BA1I? No.
    cuz ...
       i.   too slow
       ii.  mutations --> BA1B (x)
       iii. ten randomly generated DNA strings
            --> BA1I (x), solution of BA1I is only for single text
            --> concatenate ten DNA strings is indadequate (incorrect model of biological problem)

Terminology:
- motifs (k-mers)
  ; candiates for regulatory motif (or transcription factor binding site)
  ; a collection of k-mers, one from each string in DNA,
    minimizing SCORE(Motifs) among all possible choices of k-mers
- d(Pattern, text)
  ┌──────────────────────────────────────────────────────────────────────┐
  │  d(Pattern, Text) = min          HammingDistance(Pattern, Pattern')  │
  │                  all k-mers Pattern' in Text                         │
  └──────────────────────────────────────────────────────────────────────┘
- d(Pattern, Motifs)
  ┌───────────────────────────────────────────┐
  │                  t                        │
  │  d(Pattern, Dna) = Σ  d(Pattern, Dna_i)   │
  │                  i=0                      │
  └───────────────────────────────────────────┘

Plan 1.
- brute-force
- steps:
  a. candidates: all 4^k k-mer patterns as in frequency array (instead of k-mers from a string)
  b. find all (k,d)-mismatched k-mers from the first DNA string
  c. for the rest of DNA strings, do the same thing as in step (b), and intersect two sets

Plan 2.
- using BA1N neighbors function
- get len(Dnas) sets of (k,d)-mismatch k-mers, and intersect them
  TODO: improve this algorithm
        get a set of (k,d)-mismatch k-mers from the first DNA string --> candidates
        choose answers from the candidates

"""

import time


#################################################
# (k,d)-mismatched k-mers from 4^k length of array, number to pattern
#################################################

# BA01M: number to pattern
def n2p(n, k):
    """
    (int,int) -> str
    number to pattern
    """
    tbl = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    answer = ""
    while True:
        if k == 0:
            return answer
        n, k, answer = n % 4**(k - 1), k - 1, answer + tbl[n // 4**(k - 1)]


# BA01G: hamming distance
def hdist(dna1, dna2):
    """
    (str,str) -> int
    return hamming distance
    """
    return sum([dna1[i] != dna2[i] for i in range(len(dna1))])


# all at most d-mismatched k-mers
def mismatch_kmers(dna, k, d):
    """
    (str,int,int) -> {str}
    get at most d-mismatch k-mer in a DNA string
    >>> mismatch_kmers('ATT',3,1)
        {'ATC','ATG','ATA','ACT','AGT','TTT','AAT','GTT','ATT','CTT'}
    """
    answer = set()
    kmers = [n2p(i, k) for i in range(4**k)]
    for km in kmers:
        for j in range(len(dna) - k + 1):
            if hdist(dna[j:j + k], km) <= d:
                answer.add(km)
    return answer


# BA02A: motif enumeration
def motif_enum(dnas, k, d):
    """
    (str,int,int) -> {str}
    get all (k,d)-motifs in a collection of DNA strings
    TODO: How can I get the result when k and d are not given to me?
    >>> motif_enum(['ATTTGGC','TGCCTTA','CGGTATC','GAAAATT'],3,1)
        {'TTT','ATA','ATT','GTT'}
    """
    result = mismatch_kmers(dnas[0], k, d)
    for dna in dnas[1:]:
        temp = mismatch_kmers(dna, k, d)
        result = result.intersection(temp)
    return result


#################################################
# (k,d)-mismtatch k-mers from "neighbors" function
#################################################

# BA01N: neighbors
def neighbors(pattern, d):
    """
    (str,int) -> {str}
    returns all k-mers of Hamming distance at most d from Pattern.
    >>> neighbors('ACG',1)  # <-- this is immediate neighbors
        {'ACC','CCG','AGG','AAG','GCG','ATG','ACT','ACA','TCG','ACG'}
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 0:
        return {}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffixneighbors = neighbors(pattern[1:], d)
    for text in suffixneighbors:
        if hdist(text, pattern[1:]) < d:
            for nt in ['A', 'C', 'G', 'T']:
                neighborhood.add(nt + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


# BA02A: motif enumeration
def motif_enum_neighbors(dnas, k, d):
    """
    ([str],int,int) -> {str}
    get all (k,d)-motifs in a collection of DNA strings
    TODO: How can I get the result when k and d are not given to me?
    >>> motif_enum_neighbors(['ATTTGGC','TGCCTTA','CGGTATC','GAAAATT'],3,1)
        {'TTT','ATA','ATT','GTT'}
    """
    strand_kmers = [set() for _ in range(len(dnas))]
    for count,strand in enumerate(dnas):
        for i in range(len(strand) - k + 1):
            kmer = strand[i:i+k]
            strand_kmers[count].update(neighbors(kmer,d))
    result = strand_kmers[0]
    for i in range(1,len(strand_kmers)):
        result.intersection_update(strand_kmers[i])
    return result


# main
def main():
    f = open('/home/wsl/rosalind/data/ba02a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, d = [int(number.strip()) for number in lines[0].split()]
    strings = [line.strip() for line in lines[1:]]
    f.close()

    start_time = time.time()
    answer = motif_enum(strings, k, d)
    print(*answer)
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    answer = motif_enum_neighbors(strings, k, d)
    print(*answer)
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()
