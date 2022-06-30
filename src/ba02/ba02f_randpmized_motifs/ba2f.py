"""
Rosalind: BA2F
Implement Randomized Motifs Search

We will now turn to randomized algorithms that flip coins and roll dice in order to search for motifs.
Making random algorithmic decisions may sound like a disastrous idea;
just imagine a chess game in which every move would be decided by rolling a die.
However, an 18th Century French mathematician and naturalist, Comte de Buffon,
first proved that randomized algorithms are useful by randomly dropping needles onto parallel strips of wood
and using the results of this experiment to accurately approximate the constant π.

Randomized algorithms may be nonintuitive because they lack the control of traditional algorithms.
Some randomized algorithms are Las Vegas algorithms,
which deliver solutions that are guaranteed to be exact,
despite the fact that they rely on making random decisions.
Yet most randomized algorithms are Monte Carlo algorithms.
These algorithms are not guaranteed to return exact solutions,
but they do quickly find approximate solutions.
Because of their speed, they can be run many times,
allowing us to choose the best approximation from thousands of runs.

A randomized algorithm for motif finding is given below.
  ╔══════════════════════════════════════════════════════════════════════════════════════╗
  ║  RANDOMIZEDMOTIFSEARCH(Dna, k, t)                                                    ║
  ║      randomly select k-mers Motifs = (Motif1, ... , Motift) in each string from Dna  ║
  ║      BestMotifs <- Motifs                                                            ║
  ║      while forever                                                                   ║
  ║          Profile <- Profile(Motifs)         # get Profile                            ║
  ║          Motifs <- Motifs(Profile, Dna)     # find best Motifs                       ║
  ║          if Score(Motifs) < Score(BestMotifs)                                        ║
  ║              BestMotifs <- Motifs                                                    ║
  ║          else                                                                        ║
  ║              return BestMotifs                                                       ║
  ╚══════════════════════════════════════════════════════════════════════════════════════╝

Implement RandomizedMotifSearch
Given: Positive integers k and t, followed by a collection of strings Dna.

Return: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times.
        Remember to use pseudocounts!

Sample Dataset
8 5
CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
TAGTACCGAGACCGAAAGAAGTATACAGGCGT
TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
AATCCACCAGCTCCACGTGCAATGTTGGCCTA

Sample Output
TCTCGGGG
CCAAGGTG
TACAGGCG
TTCAGGTG
TCCACGTG

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
      -> GreedyMotifSearch (BA2D)
      -> PREV: GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> HERE: Randomized Motif Search (BA2F)
      -> NEXT: Gibbs Sampling (BA2G)

Chain of Functions:

  count ┌───────────────────────────────────────────────────────────── score ──┐
        └ get_profile ─── pr_score ─┬─ prof_most_probable_kmer ── get_motifs ──┼── random_motif_search
                         all_kmers ─┘                                          │
                                                               random_motifs ──┘

                     count
                 ┌─────┴────────────────┐
               get_profile            score
                 │                      │
  all_kmers    pr_score                 │
    └────────────┤                      │
         profile_most_probable_kmer     │    random_motifs
                 └──────────────────────┼──────┘
                             random_motif_search

Info:
- Review (what I have done until now)
  PROFILE(Motifs): get profile as a list of float list
  MOTIFS(Profile, Dnas): "profile-most probable k-mer"s from each DNA string

             A: 4/5   0    0   1/5        ttaccttaac
    profile  C:  0   3/5  1/5   0    DNA  gatgtctgtc
             G: 1/5  1/5  4/5   0         acggcgttag
             T:  0   1/5   0   4/5        ccctaacgag
                                          cgtcagaggt

    Motifs(Profile, Dna)   ttACCTtaac
                           gAGGTctgtc
                           acgGCGTtag
                           ccctaACGAg
                           cgtcagAGGT

- Monte-Carlo algorithm
    randomly chosen Motifs
    -> Motifs(Profile(Motifs), Dna)
    -> Profile(Motifs(Profile(Motifs), Dna))
    -> Motifs(Profile(Motifs(Profile(Motifs), Dna)))
    -> Profile(Motifs(Profile(Motifs(Profile(Motifs), Dna))))
    -> ... (a) get profile, (b) find best Motifs ... again and again

Plan 1.
- create a list of t random indices, then get random Motifs
  [idx1, idx, ... , idx_t]  -->  [DNA1[idx1:idx1+k], DNA2[idx2:idx2+k], ... , DNA_t[idx_t:idx_t + k]]

Plan 2.
- get all k-mers from each DNA strings, then randomly choose Motif from each collection of k-mers
  DNA1  --> k-mers --> Motif
  DNA2  --> k-mers --> Motif
  ...
  DNA_t --> k-mers --> Motif

═════════════════════════════════════════════════

References:
- random.randint
- random.choice
- ranodm.choices
- How do I create a list of random numbers without duplicates?
  https://stackoverflow.com/questions/9755538/how-do-i-create-a-list-of-random-numbers-without-duplicates
    import random
    random.sample(range(100), 10)   # 0..99 ten random int, without duplicate
- How do I generate a random list in Python with duplicates numbers
  https://stackoverflow.com/questions/39379515/how-do-i-generate-a-random-list-in-python-with-duplicates-numbers
    [random.choice(range(10)) for i in range(10)]
"""

import time
import random


def all_kmers(dna, k):
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


def count(motifs):
    """
    [str] -> [[int]]
    implementation of 'Count(Motifs)'
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


def score(motifs):
    """
    [str] -> int
    implementation of "Score(Motifs)"
    """
    t = len(motifs)
    return sum([t - max(ls) for ls in zip(*count(motifs))])


def get_profile(motifs, pseudocount = 1):
    """
    ([str],int) -> [[float]]
    implementation of 'Profile(Motifs)'
    """
    t = len(motifs)
    return [list(map(lambda x: (x + pseudocount) / (t*(1+pseudocount)), ls)) for ls in count(motifs)]


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
    (str,int,[[float]]) -> str
    get profile most probable k-mer from a DNA string and given profile
    ignoring the fact that there can be more than one k-mer
    >>> profile = [[0.2,0.2,0.3,0.2,0.3],[0.4,0.3,0.1,0.5,0.1],[0.3,0.3,0.5,0.2,0.4],[0.1,0.2,0.1,0.1,0.2]]
    >>> dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
    >>> prof_most_probable_kmer(dna,5,profile)
        'CCGAG'   # <-- modified to get a string, not set of strings
    """
    result = {}
    largest = 0
    for kmer in all_kmers(dna, k):
        val = pr_score(profile, kmer)
        result[kmer] = val
        if val > largest:
            largest = val
    for k, v in result.items():
        if v == largest:
            return k


def get_motifs(dnas, k, profile):
    """
    ([[float]],[str]) -> [str]
    implementation of 'MOTIFS(Profile, Dna)' function
    returns a collection of profile-most probable k-mers from each DNA string
    >>> profile = [[0.2,0.2,0.3,0.2,0.3],[0.4,0.3,0.1,0.5,0.1],[0.3,0.3,0.5,0.2,0.4],[0.1,0.2,0.1,0.1,0.2]]
    >>> dnas = ['CGCCCC','GGGCGA','TAGTAC','TAGATC','AATCCA']
    >>> get_motifs(dnas,5,profile)
        ['CGCCC','GGGCG','TAGTA','TAGAT','ATCCA']
    """
    motifs = []
    for dna in dnas:
        motif = prof_most_probable_kmer(dna, k, profile)
        motifs.append(motif)
    return motifs


def random_motifs_all_kmers(dnas,k):
    """
    ([str],int) -> [str]
    returns a collection of randomly selected k-mer Motifs from DNA string
    all k-mers from each DNA string, then select each Motif
    Getting all kmers from each DNA string is time-consuming process. (TODO: improve this)
    >>> dnas = ['CGCCCC','GGGCGA','TAGTAC','TAGATC','AATCCA']
    >>> random_motifs_all_kmers(dnas,3)
        ['GCC','CGA','TAG','AGA','AAT']
    """
    rand_motifs = []
    for dna in dnas:
        kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]
        motif = random.choice(kmers)
        rand_motifs.append(motif)
    return rand_motifs


def random_motifs(dnas, k):
    """
    ([str],int) -> [str]
    returns a collection of randomly selected k-mer Motifs from each DNA string
    just random indices, then random k-mers
    >>> dnas = ['CGCCCC','GGGCGA','TAGTAC','TAGATC','AATCCA']
    >>> random_motifs(dnas,3)
        ['CCC','CGA','AGT','GAT','TCC']  # <-- not fixed result
    """
    n = len(dnas[0])
    t = len(dnas)
    rand_indices = [random.choice(range(n - k + 1)) for i in range(t)]
    rand_motifs = []
    for i in range(t):
        motif = dnas[i][rand_indices[i]:rand_indices[i]+k]
        rand_motifs.append(motif)
    return rand_motifs


def random_motif_search(dnas, k, t):
    """
    ([str],int,int) -> [str]
    implementation of 'RANDOMIZEDMOTIFSEARCH' function
    >>> dnas = ['CGCCCC','GGGCGA','TAGTAC','TAGATC','AATCCA']
    >>> random_motif_search(dnas,3,5)
        ['CGC','GGC','AGT','ATC','ATC']
    """
    best_motifs = random_motifs(dnas, k)
    while True:
        profile = get_profile(best_motifs)
        motifs = get_motifs(dnas, k, profile)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def repeat(DNAs, k, niter):
    """
    ([str],int,int) -> [str]
    """
    best_score = float('inf')
    best_motifs = []
    t = len(DNAs)
    for i in range(niter):
        motifs = random_motif_search(DNAs, k, t)
        cur_score = score(motifs)
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
    return best_motifs


def main():
    f = open('/home/wsl/rosalind/data/ba02f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, t = [int(num.strip()) for num in lines[0].split()]
    dnas = [line.strip() for line in lines[1:]]
    f.close()


    start_time = time.time()
    answer = repeat(dnas,k,1000)
    for motif in answer:
        print(motif)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
