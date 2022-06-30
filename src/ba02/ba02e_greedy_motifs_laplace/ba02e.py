"""
Rosalind: BA2E
Implement GreedyMotifSearch (pseudocount)

We encountered GreedyMotifSearch in “Implement GreedyMotifSearch”.
In this problem, we will power it up with pseudocounts.

Implement GreedyMotifSearch with Pseudocounts
Given: Integers k and t, followed by a collection of strings Dna.

Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t) with pseudocounts.
If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

Sample Dataset
3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG

Sample Output
TTC
ATC
TTC
ATC
TTC

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
      -> Median string problem (BA2B) - not fast enough... find another algorithm
      ↓
    * GreedyMotifSearch
      -> Profile-most probable k-mer problem (BA2C)
      -> PREV: GreedyMotifSearch (BA2D)
      -> HERE: GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> NEXT: Randomized Motif Search (BA2F)

Info:
- Why do I need pseudo-count (Laplace’s Rule of Succession)?

                          A:  .2   .2  .0   .0   .0   .0   .9   .1   .1   .1   .3   .0
                          C:  .1   .6  .0   .0   .0   .0   .0   .4   .1   .2   .4   .6
                          G:  .0   .0   1    1   .9   .9   .1   .0   .0   .0   .0   .0
                          T:  .7   .2  .0   .0   .1   .1   .0   .5   .8   .7   .3   .4
  Pr(TCGTGGATTTCC|Profile) =  .7 * .6 * 1 * .0 * .9 * .9 * .9 * .5 * .8 * .7 * .4 * .6 = 0 (ZERO!!)

- applying Laplace’s Rule of Succession (Pseudo-count)
  Motifs            T    A    A    C
                    G    T    C    T
                    A    C    T    A
                    A    G    G    T

  COUNT(Motifs) A:  2    1    1    1       2+1   1+1   1+1   1+1
                C:  0    1    1    1  ->   0+1   1+1   1+1   1+1
                G:  1    1    1    0       1+1   1+1   1+1   0+1
                T:  1    1    1    2       1+1   1+1   1+1   2+1

  PROFILE(Motifs) 2/4  1/4  1/4  1/4       3/8   2/8   2/8   2/8
                    0  1/4  1/4  1/4  ->   1/8   2/8   2/8   2/8
                  1/4  1/4  1/4    0       2/8   2/8   2/8   1/8
                  1/4  1/4  1/4  2/4       2/8   2/8   2/8   3/8

Plan 1.
"""

import time


def count(motifs):
    """
    [str] -> [[int]]
    Count(Motifs)
    t * k matrix of nucleotide counting (row-by-row)
    """
    dict = {'A':[], 'C':[], 'G':[], 'T':[]}
    for k in range(len(motifs[0])):
        temp = {'A':0, 'C':0, 'G':0, 'T':0}
        for t in range(len(motifs)):
            temp[motifs[t][k]] += 1
        for k,v in temp.items():
            dict[k].append(v)
    result = []
    for nuc in "ACGT":
        result.append(dict[nuc])
    return result


def get_profile(motifs, pseudocount = 1):
    """
    ([str],int) -> [[float]]
    Profile(Motifs)
    t * k matrix of nucleotide probability (row-by-row)
    """
    t = len(motifs)
    return [list(map(lambda x: (x + pseudocount) / (t*(1+pseudocount)), ls)) for ls in count(motifs)]


def score(motifs):
    """
    [str] -> int
    Score(Motifs)
    """
    t = len(motifs)
    return sum([t - max(ls) for ls in zip(*count(motifs))])


def pr_score(profile, pattern):
    """
    ([[float]],str) -> float
    Pr(Pattern|Profile)
    """
    val = 1
    for i in range(len(pattern)):
        nuc = pattern[i]
        if nuc =='A':
            val *= profile[0][i]
        elif nuc == 'C':
            val *= profile[1][i]
        elif nuc == 'G':
            val *= profile[2][i]
        elif nuc =='T':
            val *= profile[3][i]
    return val


def all_kmers(dna, k):
    """
    (str,int) -> [str]
    get all k-mers from a string
    TODO: What if I return generator instead of list?
    """
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


def all_kmers_generator(dna, k):
    """
    (str,int) -> gen
    get all k-mers from a string - generator version
    """
    return (dna[i:i+k] for i in range(len(dna) - k + 1))


def greedy_motif_search_pseudocount(dna, k, t):
    """
    (str,int,int) -> [str]
    version 1. greedy motifs search with pseudocount
    """
    best_motifs = [s[:k] for s in dna]
    n = len(dna[0])
    for i in range(n - k + 1):
        motif = dna[0][i:i+k]
        motifs = [motif]
        for j in range(1, t):
            pros = get_profile(motifs)
            max_prob_score = 0
            max_k_mer = dna[j][:k]
            for k_mer in all_kmers(dna[j], k):
                k_mer_score = pr_score(pros, k_mer)
                if k_mer_score > max_prob_score:
                    max_prob_score = k_mer_score
                    max_k_mer = k_mer
            motifs.append(max_k_mer)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def greedy_motif_search_pseudocount_generator(dna, k, t):
    """
    (str,int,int) -> [str]
    version 2. greedy motifs search with pseudocount
    """
    best_motifs = [s[:k] for s in dna]
    n = len(dna[0])
    for i in range(n - k + 1):
        motif = dna[0][i:i+k]
        motifs = [motif]
        for j in range(1, t):
            pros = get_profile(motifs)
            max_prob_score = 0
            max_k_mer = dna[j][:k]
            for k_mer in all_kmers_generator(dna[j], k):
                k_mer_score = pr_score(pros, k_mer)
                if k_mer_score > max_prob_score:
                    max_prob_score = k_mer_score
                    max_k_mer = k_mer
            motifs.append(max_k_mer)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def main():
    f = open('/home/wsl/rosalind/data/ba02e.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, t = [int(num.strip()) for num in lines[0].split()]
    dnas = [line.strip() for line in lines[1:]]
    f.close()

    start_time = time.time()
    for s in greedy_motif_search_pseudocount(dnas, k, t):
        print(s)
    print("--- %s seconds ---" % (time.time() - start_time))  # --- 1.3414032459259033 seconds ---

    start_time = time.time()
    for s in greedy_motif_search_pseudocount_generator(dnas, k, t):
        print(s)
    print("--- %s seconds ---" % (time.time() - start_time))  # --- 1.3130505084991455 seconds ---


if __name__ == "__main__":
    main()