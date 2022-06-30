"""
Rosalind: BA2D
Implement GreedyMotifSearch

  ╔═════════════════════════════════════════════════════════════════════════════════╗
  ║  GREEDYMOTIFSEARCH(Dna, k, t)                                                   ║
  ║      BestMotifs <- motif matrix formed by first k-mers in each string from Dna  ║
  ║      for each k-mer Motif in the first string from Dna                          ║
  ║          Motif1 <- Motif                                                        ║
  ║          for i = 2 to t                                                         ║
  ║              form Profile from motifs Motif1, ..., Motif_i-1                    ║
  ║              Motifi <- Profile-most probable k-mer in the i-th string in Dna    ║
  ║          Motifs <- (Motif1, ..., Motift)                                        ║
  ║          if SCORE(Motifs) < SCORE(BestMotifs)                                   ║
  ║              BestMotifs <- Motifs                                               ║
  ║      return BestMotifs                                                          ║
  ╚═════════════════════════════════════════════════════════════════════════════════╝

Implement GreedyMotifSearch
Given: Integers k and t, followed by a collection of strings Dna.

Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t).
If at any step you find more than one Profile-most probable k-mer in a given string,
use the one occurring first.

Sample Dataset
3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG

Sample Output
CAG
CAG
CAA
CAA
CAA

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
      -> PREV: Profile-most probable k-mer problem (BA2C)
      -> HERE: GreedyMotifSearch (BA2D)
      -> NEXT: GreedyMotifSearch with pseudo-count (BA2E)

Plan 1.
- implementing GreedyMotifSearch using "Profile-most probable k-mer" function in BA2C
- Pseudocode:
  ╔═════════════════════════════════════════════════════════════════════════════════╗
  ║  GREEDYMOTIFSEARCH(Dna, k, t)                                                   ║
  ║      BestMotifs <- motif matrix formed by first k-mers in each string from Dna  ║
  ║      for each k-mer Motif in the first string from Dna                          ║
  ║          Motif1 <- Motif                                                        ║
  ║          for i = 2 to t                                                         ║
  ║              form Profile from motifs Motif1, ..., Motif_i-1                    ║
  ║              Motifi <- Profile-most probable k-mer in the i-th string in Dna    ║
  ║          Motifs <- (Motif1, ..., Motift)                                        ║
  ║          if SCORE(Motifs) < SCORE(BestMotifs)                                   ║
  ║              BestMotifs <- Motifs                                               ║
  ║      return BestMotifs                                                          ║
  ╚═════════════════════════════════════════════════════════════════════════════════╝

- practice with sample dataset
  not easy, divide algorithm into several steps
  (1) Let's make two assumptions.
      One for "best Motifs", the other for best "Motif"
      Motifs = {Motif1, Motif2, ... , Motif_t}

      Assumption 1. The best Motifs == first k-mers in each DNA strings (len(Motifs) == t)
      Assumption 2. The best Motif in a DNA string == first k-mer

                              ↓first                              ↓last
      GGCGTTCAGGCA --> kmers: GGC GCG CGT GTT TTC TCA CAG AGG GGC GCA
      AAGAATCAGTCA            one of these k-mers will be a Motif in the first DNA string
      CAAGGAGTTCGC            Let's suppose the Motif is the first k-mer 'GGC' (by Assumption 1.),
      CACGTCAATCAC            and best Motifs is ['GGC','AAG','CAA','CAC','CAA'].
      CAATAATATTCG

  (2) Put the best Motifs (Motifs_0: first t k-mers (t * Motif_0)) aside for a while.

  (3) Find the best Motifs (Motifs_i) for each i-th k-mer (Motif_i) in the first DNA string.
      And get the NEW best Motifs by comparing SCORE(Motifs_0) with SCORE(Motifs_i).
      Do this process again and again till the last k-mer in the first DNA string.
      I can get the real BEST Motifs by this "len(DNA[0]) - k + 1" iteration.

═════════════════════════════════════════════════

References:
- StackOverFlow: Transpose list of lists
  https://stackoverflow.com/questions/6473679/transpose-list-of-lists
- StackOverFlow: Get the cartesian product of a series of lists?
  https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists
- use generator if possible!
- Greedy Algorithm
  https://www.programiz.com/dsa/greedy-algorithm#:~:text=A%20greedy%20algorithm%20is%20an,if%20the%20choice%20is%20wrong.
"""

import time
import math


#################################################
# Just Practice: Finding Consensus string from given DNA strings
#################################################


def count_cbc(motifs):
    """
    ([str]) -> [[int]]
    implementation of "COUNT(Motifs)" function
    list of list, count column-by-column
    >>> motifs = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    >>> count_cbc(motifs)
        [[1,3,1,0],[4,0,1,0],[2,2,1,0],[1,0,3,1],[2,0,1,2],[2,1,0,2],[1,2,1,1],[4,0,0,1],[0,0,2,3],[0,2,1,2],[1,3,1,0],[2,2,1,0]]
    """
    lls = []  # list of list of four integer values
    for tup in zip(*motifs):
        ls = [0,0, 0, 0]  # A, C, G, T
        for nuc in tup:
            if nuc == 'A':
                ls[0] += 1
            elif nuc == 'C':
                ls[1] += 1
            elif nuc == 'G':
                ls[2] += 1
            elif nuc == 'T':
                ls[3] += 1
        lls.append(ls)
    return lls


def score_cbc(motifs):
    """
    ([str]) -> int
    implementation of "SCORE(Motifs): function
    >>> motifs = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    >>> score_cbc(motifs)
        28
    """
    return sum([len(motifs) - max(ls) for ls in count_cbc(motifs)])


def profile_cbc(motifs):
    """
    ([str]) -> [[float]]
    profile(Motifs), column-by-column
    It is better to create a list with the same column of multiple strings.
    >>> motifs = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    >>> profile_cbc(motifs)
        [[0.2,0.6,0.2,0.0],[0.8,0.0,0.2,0.0],[0.4,0.4,0.2,0.0],[0.2,0.0,0.6,0.2],[0.4,0.0,0.2,0.4],[0.4,0.2,0.0,0.4],[0.2,0.4,0.2,0.2],[0.8,0.0,0.0,0.2],[0.0,0.0,0.4,0.6],[0.0,0.4,0.2,0.4],[0.2,0.6,0.2,0.0],[0.4,0.4,0.2,0.0]]
    """
    t = len(motifs)
    return [list(map(lambda x: x / t,ls)) for ls in count_cbc(motifs)]


def profile_entropy_cbc(motifs):
    """
    ([str]) -> [[float]]
    profile of entropy
    -(pr1 * log2(pr1) + pr2*log2(pr2) + ... + prN*log2(prN))
    >>> motifs = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    >>> profile_entropy_cbc(motifs)
        [1.3709505944546687, 0.7219280948873623, 1.5219280948873621, 1.3709505944546687, 1.5219280948873621, 1.5219280948873621, 1.9219280948873623, 0.7219280948873623, 0.9709505944546686, 1.5219280948873621, 1.3709505944546687, 1.5219280948873621]
    """
    ls = []
    for lpr in profile_cbc(motifs):
        val = 0
        for pr in lpr:
            if pr == 0:
                val += 0
            else:
                val += pr * math.log2(pr)
        ls.append(-val)
    return ls


def product(ar_list):
    """
    ([[a]]) -> [[a]]
    returns cartesian product from list of list
    >>> list(product([[1,2,3],[4,5,6]]))
        [[1,4],[1,5],[1,6],[2,4],[2,5],[2,6],[3,4],[3,5],[3,6]]
    """
    if not ar_list:
        yield []
    else:
        for a in ar_list[0]:
            for prod in product(ar_list[1:]):
                yield [a,] + prod


def consensus_cbc(motifs):
    """
    ([str]) -> [str]
    implementation of "CONSENSUS(Motifs)" function
    Consider that there can be more than one,
    because there can be multiple values with the same max-frequency value.
    So we need 'product' function.
    """
    # find index of the largest probability in each profile matrix
    # consider that there can be more than one!!
    ls_freq = []
    for ls in profile_cbc(motifs):
        max_pr = max(ls)
        temp = []
        for i in range(len(ls)):  # if two or more nucleotides have max pr-value
            if ls[i] == max_pr:
                temp.append(i)
        ls_freq.append(temp)
    # get consensus motifs from above result
    prod = product(ls_freq)
    candidates = []
    for ls in prod:  # Be careful not to use list itself. Use generator here!
        temp_str = ""
        for num in ls:
            if num == 0:
                temp_str += 'A'
            elif num == 1:
                temp_str += 'C'
            elif num == 2:
                temp_str += 'G'
            elif num == 3:
                temp_str += 'T'
        candidates.append(temp_str)
    return candidates


def pr_score_rbr(profile, kmer):
    """
    ([[float],str]) -> float
    returns probability score by row-by-row
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


#################################################
# Finding best motifs (greedy motif searching)
#################################################


def count(motifs):
    """
    [str] -> {char:[int]}
    implementation of "Count(Motifs)" function
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


def get_profile(motifs):
    """
    [str] -> [[float]]
    implementation of "Profile(Motifs)"
    t * k matrix of nucleotide probability (row-by-row)
    TODO: Find the problem of this function!
          If count value is zero, what will happen?
    """
    t = len(motifs)
    return [list(map(lambda x: x / t, ls)) for ls in count(motifs)]


def score(motifs):
    """
    [str] -> int
    implementation of "Score(Motifs)"
    """
    t = len(motifs)
    return sum([t - max(ls) for ls in zip(*count(motifs))])


def pr_score(profile, pattern):
    """
    ([[float]],str) -> float
    implementation of "Pr(Pattern|Profile)"
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
    """
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


# BA02D: greed motif search
def greedy_motif_search(dnas, k, t):
    """
    (str,int,int) -> [str]
    greedy motif search, algorithm in textbook
    >>> dnas = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    >>> greedy_motif_search(dnas,3,5)
        ['CAG','CAG','CAA','CAA','CAA']
    """
    # Let's suppose that a collection of the first k-mers from each DNA string is the "Motifs".
    best_motifs = [s[:k] for s in dnas]
    n = len(dnas[0])
    # One of the k-mer from the first DNA string must be the first "Motif" (That's true!) --> loop!
    for i in range(n - k + 1):
        motif = dnas[0][i:i+k]
        motifs = [motif]
        # Find the "Motifs" when the first k-mer "Motif" is DNA[0][i:i+k].
        for j in range(1, t):
            pros = get_profile(motifs)
            max_prob_score = 0
            max_k_mer = dnas[j][:k]
            # find a k-mer Motif with max probability score from j-th DNA string
            for k_mer in all_kmers(dnas[j], k):
                k_mer_score = pr_score(pros, k_mer)
                if k_mer_score > max_prob_score:
                    max_prob_score = k_mer_score
                    max_k_mer = k_mer
            motifs.append(max_k_mer)
        # find new best Motifs
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def main():
    f = open('/home/wsl/rosalind/data/ba02d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, t = [int(num.strip()) for num in lines[0].split()]
    dnas = [line.strip() for line in lines[1:]]
    f.close()

    start_time = time.time()
    for motif in greedy_motif_search(dnas, k, t):
        print(motif)  # --- 1.3052973747253418 seconds ---
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
