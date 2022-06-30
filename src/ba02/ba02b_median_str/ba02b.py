"""
Rosalind: BA2B
Find a Median String

Given a k-mer Pattern and a longer string Text,
we use d(Pattern, Text) to denote the minimum Hamming distance
between Pattern and any k-mer in Text,

d(Pattern, Text) = min                 HammingDistance(Pattern, Pattern')
                  all k-mers Pattern' in Text

    example:
    d(GATTCTCA, gcaaaGACGCTGAccaa) = 3

Given a k-mer Pattern and a set of strings
Dna = {Dna_1, ... , Dna_t},
we define d(Pattern, Dna) as the sum of distances
between Pattern and all strings in Dna,

                    t
    d(Pattern, Dna) = Σ  d(Pattern, Dna_i)
                    i=0

    example:
           d(AAA, Dna) = 1 + 1 + 2 + 0 + 1 = 5
           ttaccttAAC   1
           gATAtctgtc   1
    Dna    ACGgcgttcg   2
           ccctAAAgag   0
           cgtcAGAggt   1

Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna) over all
k-mers Pattern, the same task that the Equivalent Motif Finding Problem is
trying to achieve. We call such a k-mer a median string for Dna.

Median String Problem
Find a median string.

Given: An integer k and a collection of strings Dna.

Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern.
(If multiple answers exist, you may return any one.)

Sample Dataset
3
AAATTGACGCAT
GACGACCACGTT
CGTCAGCGCCTG
GCTGAGCACCGG
AGTACGGGACAG

Sample Output
GAC

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> PREV: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference? (very important)
      -> HERE: "Median string problem" (BA2B) - not fast enough... find better algorithm
      ↓
    * GreedyMotifSearch
      -> NEXT: Profile-most probable k-mer problem (BA2C)

Info:
- Counting row-by-row returns exactly the same result as counting col-by-col.
                                                               (2)
                          T  C  G  G  G  G  g  T  T  T  t  t   3 row-by-row
                          c  C  G  G  t  G  A  c  T  T  a  C   4
                          a  C  G  G  G  G  A  T  T  T  t  C   2
                          T  t  G  G  G  G  A  c  T  T  t  t   4
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C   3
                          T  t  G  G  G  G  A  c  T  T  C  C   2
                          T  C  G  G  G  G  A  T  T  c  a  t   3
                          T  C  G  G  G  G  A  T  T  c  C  t   2
                          T  a  G  G  G  G  A  a  c  T  a  C   4
                          T  C  G  G  G  t  A  T  a  a  C  C + 3   Interesting!
    SCORE(Motifs)   (1)   3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30  result score is the same!
                          col-by-col

   (1) given motifs -> consensus
      ; Motifs -> CONSENSUS(Motifs)
      ; Motif Finding Problem
   (2) condidates of consensus -> best possible collection Motifs for each consensus string
      ; CONSENSUS(Motifs) -> Motifs
      ; Equivalent Motif Finding Problem

                      (1) first approach             (2) second approach
                      T C G G G G g T T T t t        T C G G G G g T T T t t   3
                      c C G G t G A c T T a C        c C G G t G A c T T a C   4
                      a C G G G G A T T T t C        a C G G G G A T T T t C   2
                      T t G G G G A c T T t t        T t G G G G A c T T t t   4
    Motifs            a a G G G G A c T T C C        a a G G G G A c T T C C   3
                      T t G G G G A c T T C C        T t G G G G A c T T C C   2
                      T C G G G G A T T c a t        T C G G G G A T T c a t   3
                      T C G G G G A T T c C t        T C G G G G A T T c C t   2
                      T a G G G G A a c T a C        T a G G G G A a c T a C   4
                      T C G G G t A T a a C C        T C G G G t A T a a C C + 3
    SCORE(Motifs)     3+4+0+0+1+1+1+5+2+3+6+4 = 30   3+4+0+0+1+1+1+5+2+3+6+4= 30

- comparing three problems
  ┌──────────────────────────────────────────────────────────────────┐
  │ Motif Finding Problem:                                           │ O(n^t * k * t)
  │   Given a collection of strings,                                 │
  │   find a set of k-mers, one from each string,                    │
  │   that minimizes the score of the resulting motif.               │
  │                                                                  │
  │ Input:                                                           │
  │   A collection of strings Dna and an integer k.                  │
  │ Output:                                                          │
  │   A collection "Motifs" of k-mers, one from each string in Dna,  │ Motifs
  │   minimizing SCORE(Motifs) among all possible choices of k-mers. │
  └──────────────────────────────────────────────────────────────────┘
  ┌────────────────────────────────────────────────────────────────────┐
  │ Equivalent Motif Finding Problem:                                  │
  │   Given a collection of strings,                                   │
  │   find a pattern and a collection of k-mers (one from each string) │
  │   that minimizes the distance between all possible patterns        │
  │   and all possible collections of k-mers.                          │
  │                                                                    │
  │ Input:                                                             │
  │   A collection of strings Dna and an integer k.                    │
  │ Output:                                                            │
  │   A k-mer "Pattern" and a collection of k-mers "Motifs",           │ Pattern, Motifs
  │   one from each string in Dna, minimizing d(Pattern, Motifs)       │
  │   among all possible choices of Pattern and Motifs.                │
  └────────────────────────────────────────────────────────────────────┘
  ┌──────────────────────────────────────────────────┐
  │ Median String Problem:                           │
  │   Find a median string.                          │
  │ Input:                                           │
  │   A collection of strings Dna and an integer k.  │
  │ Output:                                          │
  │   A k-mer Pattern minimizing d(Pattern, Dna)     │
  │   among all k-mers Pattern.                      │
  └──────────────────────────────────────────────────┘

- What is the meaning of above comparison?
  ; The col-by-col approach is exactly the same as "Motif finding problem", nothing else.
  ; But row-by-row approach is something different.
    The key observation for solving the Equivalent Motif Finding Problem is that,
    given Pattern, we don’t need to explore all possible collections (4^k) Motifs
    in order to minimize d(Pattern, Motifs).

    MOTIFS(AAA, Dna):
                         ttaccttAAC
                         gATAtctgtc
                     Dna ACGgcgttcg
                         ccctAAAgag
                         cgtcAGAggt

Plan 1.
- What is the "median string"?
  ; a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern

- What are the meanings of "d(pattern, text)" and "d(pattern, Dna)"?
  d(pattern, one Dna string)       = x   (minimum hamming distance value)
  d(pattern, multiple Dna strings) = x + y + ... + z

-

Plan 2.
- pseudocode:
  ╔═══════════════════════════════════════════════════════╗
  ║  MEDIANSTRING(Dna, k)                                 ║ O(4^k * n * k * t)
  ║      distance <- infinite                             ║
  ║      for each k-mer Pattern from AA...AA to TT...TT   ║
  ║          if distance > d(Pattern, Dna)                ║
  ║              distance <- d(Pattern, Dna)              ║
  ║              Median <- Pattern                        ║
  ║      return Median                                    ║
  ╚═══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

"""

import time


def n2p(n, k):
    tbl = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    answer = ""
    while True:
        if k == 0:
            return answer
        n, k, answer = n % 4**(k - 1), k - 1, answer + tbl[n // 4**(k - 1)]
    return answer


def hdist(dna1, dna2):
    """
    returns int - hamming distance
    """
    return sum([dna1[i] != dna2[i] for i in range(len(dna1))])


def all_kmers(text, k):
    """
    returns [string] - all k-mers from text
    """
    return [text[i:i+k] for i in range(len(text)-k+1)]


def min_hdist(pattern, kmers):
    """
    find minimum hamming distance with a pattern and all kmers
    """
    return min([hdist(pattern, kmer) for kmer in kmers])


# all possible k-mer patterns (4^k)
def all_patterns(k):
    """
    int -> [str]
    get all k-characters patterns - 4^k
    """
    return [n2p(i, k) for i in range(4**k)]


# BA02B: median strings
def median_string(lstr, k):
    """
    (int,[str]) -> {str}
    """
    min_val = float('inf')
    patterns = set()
    for pat in all_patterns(k):
        sum_val = 0
        for str in lstr:
            sum_val += min_hdist(pat, all_kmers(str, k))
        if sum_val < min_val:
            min_val = sum_val
            patterns = {pat}
        elif sum_val == min_val:
            patterns.add(pat)
    return patterns


#################################################
# using textbook algorithm as similar as possible
#################################################

# 'd(Pattern,Dna)' function
def d_pattern_dnas(pattern, dnas):
    """
    (str,[str]) -> {str}
    implementing 'd(Pattern,Dna)' function
    >>> d_pattern_dnas("AAA",['TTACCTTAAC','GATATCTGTC','ACGGCGTTCG','CCCTAAAGAG','CGTCAGAGGT'])
        5
    """
    k = len(pattern)
    result = 0
    for dna in dnas:
        min_val = float('inf')
        kmers = all_kmers(dna,k)
        for kmer in kmers:
            distance = hdist(pattern, kmer)
            if min_val > distance:
                min_val = distance
        result += min_val
    return result


# BA02B: median strings
def median_string_textbook(dnas, k):
    """
    ([str],int) -> {str}
    >>> median_string_textbook(['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTACGGGACAG'],3)
        {'ACG'}
    """
    distance = float('inf')
    median = set()
    for pattern in [n2p(i, k) for i in range(4**k)]:
        dist = d_pattern_dnas(pattern, dnas)
        if distance > dist:
            distance = dist
            median = {pattern}
        # elif distance == dist:      # <-- if you want all median strings, uncomment these two lines
        #     median.update({pattern})
    return median


# main
def main():
    f = open('/home/wsl/rosalind/data/ba02b.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0].strip())
    dnas = lines[1:]
    f.close()

    start_time = time.time()
    answers = median_string_textbook(dnas, k)
    for answer in answers:
        print(answer)
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    answers = median_string(dnas, k)
    for answer in answers:
        print(answer)
    print("--- %s seconds ---" % (time.time() - start_time)) # --- 2.696305751800537 seconds ---


# execute
if __name__ == "__main__":
    main()
