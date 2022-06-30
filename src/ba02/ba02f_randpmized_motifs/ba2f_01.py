"""
Rosalind: BA2F
Implement GreedyMotifSearch

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

"""

from cmath import inf
import time
import random


def count(motifs):
    """
    [str] -> [[int]]
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
    """
    t = len(motifs)
    return [list(map(lambda x: (x + pseudocount) / (t*(1+pseudocount)), ls)) for ls in count(motifs)]


def score(motifs):
    """
    [str] -> int
    """
    t = len(motifs)
    return sum([t - max(ls) for ls in zip(*count(motifs))])


def pr_score(profile, pattern):
    """
    ([[float]],str) -> float
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
    """
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


def get_motifs(profile, DNAs):
    """
    ([[float]],[str]) -> [str]
    Motifs(profile, DNAs)
    """
    motifs = []
    k = len(profile[0])
    for string in DNAs:
        kmers = all_kmers(string, k)
        max_prob = 0
        best_kmer = kmers[0]
        for kmer in kmers:
            prob = pr_score(profile, kmer)
            if prob > max_prob:
                max_prob = prob
                best_kmer = kmer
        motifs.append(best_kmer)
    return motifs


def rand_motif_search_wrong(DNAs, k, t):
    """
    ([str],int,int) -> [str]
    """
    # TODO: Why is this wrong?
    # create a randomly chosen Motifs
    best_motifs = []
    for i in range(t):
        idx_start = random.randint(0, len(DNAs[0]) - k)  # <-- This part is wrong
        randomly_selected_motif = DNAs[i][idx_start:idx_start+k]
        best_motifs.append(randomly_selected_motif)
    # find best motifs
    while True:
        profile = get_profile(best_motifs)
        new_motifs = get_motifs(profile, best_motifs)
        if score(new_motifs) < score(best_motifs):
            best_motifs = new_motifs
        else:
            return best_motifs


def rand_motif_search_wrong2(DNAs, k, t):
    """
    ([str],int,int) -> [str]
    """
    # TODO: Why is this wrong?
    rand_indices = []
    for i in range(t):
        idx_start = random.randint(0, len(DNAs[0]) - k)
        rand_indices.append(idx_start)
    rand_motifs = []
    for i in range(t):
        start = rand_indices[i]
        rand_motifs.append(DNAs[i][start:start+k])
    best_motifs = rand_motifs
    while True:
        profile = get_profile(best_motifs)
        new_motifs = get_motifs(profile, best_motifs)
        if score(new_motifs) < score(best_motifs):
            best_motifs = new_motifs
        else:
            return best_motifs


def rand_motif_search(DNAs, k, t):
    """
    ([str],int,int) -> [str]
    """
    DNAs_k_mers = [all_kmers(dna, k) for dna in DNAs]
    # dna_k_mers = map(all_kmers, dna, [k] * len(dna))
    motifs = [random.choice(k_mers) for k_mers in DNAs_k_mers]
    best_motifs = motifs
    while True:
        profile = get_profile(motifs)
        motifs = get_motifs(profile, DNAs)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def repeat(DNAs, k, niter):
    """
    ([str],int,int) -> [str]
    """
    best_score = inf
    best_motifs = []
    t = len(DNAs)
    for i in range(niter):
        motifs = rand_motif_search(DNAs, k, t)
        cur_score = score(motifs)
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
    return best_motifs


def main():
    f = open('/home/wsl/rosalind/data/ba02f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, t = [int(num.strip()) for num in lines[0].split()]
    DNAs = [line.strip() for line in lines[1:]]
    f.close()


    answer = rand_motif_search_wrong(DNAs,k,t)
    print(answer)


    start_time = time.time()
    for i in repeat(DNAs, k, 10):
        print(i)  # --- 48.222848892211914 seconds ---
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
