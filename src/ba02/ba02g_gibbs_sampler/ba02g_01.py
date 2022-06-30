"""
Rosalind: BA2G
Implement GibbsSampler

We have previously defined the notion of a Profile-most probable k-mer in a string.
We now define a Profile-randomly generated k-mer in a string Text.
For each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile),
resulting in n = |Text| - k + 1 probabilities (p1, …, pn).
These probabilities do not necessarily sum to 1,
but we can still form the random number generator Random(p1, …, pn) based on them.
GIBBSSAMPLER uses this random number generator
to select a Profile-randomly generated k-mer at each step:
if the die rolls the number i,
then we define the Profile-randomly generated k-mer as the i-th k-mer in Text.

Pseudocode:
  ╔═══════════════════════════════════════════════════════════════════════════════════╗
  ║ GIBBSSAMPLER(Dna, k, t, N)                                                        ║
  ║     randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna ║
  ║     BestMotifs <- Motifs                                                          ║
  ║     for j <- 1 to N                                                               ║
  ║         i <- Random(t)                                                            ║
  ║         Profile <- profile matrix constructed from all strings in Motifs          ║
  ║                    except for Motifi                                              ║
  ║         Motifi <- Profile-randomly generated k-mer in the i-th sequence           ║
  ║         if Score(Motifs) < Score(BestMotifs)                                      ║
  ║             BestMotifs <- Motifs                                                  ║
  ║     return BestMotifs                                                             ║
  ╚═══════════════════════════════════════════════════════════════════════════════════╝

Implement GibbsSampler
Given: Integers k, t, and N, followed by a collection of strings Dna.

Return: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts.
Remember to use pseudocounts!

Sample Dataset
8 5 100
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
import random
import time


def all_kmers(dna, k):
    """
    (str,int) -> [str]
    get all k-mers from a string
    """
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]


def count(motifs):
    """
    [str] -> [[int]]
    Count(Motifs)
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
    Score(Motifs)
    """
    t = len(motifs)
    return sum([t - max(ls) for ls in zip(*count(motifs))])


def get_profile(motifs, pseudocount = 1):
    """
    ([str],int) -> [[float]]
    """
    t = len(motifs[0])
    return [list(map(lambda x: (x+pseudocount) / (t+t*pseudocount), ls)) for ls in count(motifs)]


def pr_score(profile, pattern):
    """
    ([[float]],str) -> float
    Pr(pattern | profile)
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


def lst_probs(profile, kmers):
    """
    ([[float]],[str]) -> [float]
    get a list of probabilities of all k-mers from a string
    """
    probs = [pr_score(profile, kmer) for kmer in kmers]
    total = sum(probs)
    return [pr / total for pr in probs]


def random_idx(lst_probs):
    """
    [float] -> float
    Random(p1, p2, ..., pn)
    """
    return random.choices(range(len(lst_probs), lst_probs, k=1))[0]


def profile_randomly_gen_kmer_verbose(profile, dna):
    """
    ([[float]],str) -> str
    profile randomly generated k-mer
    """
    k = len(profile[0])
    kmers = all_kmers(dna, k)
    probs = lst_probs(profile, kmers)
    idx = random_idx(probs)
    return kmers[idx]


def profile_randomly_gen_kmer(profile, dna):
    """
    ([[float]],str) -> str
    profile randomly generated k-mer, compact version
    """
    k = len(profile[0])
    kmers = all_kmers(dna, k)
    n = len(kmers)
    probs = list(map(pr_score, [profile] * n, kmers))
    c = sum(probs)
    norm_probs = [prob / c for prob in probs]
    motif = random.choices(kmers, weights=norm_probs, k=1)
    return motif.pop()


def gibbs_sampler(dnas, k, t, n):
    """
    ([str],int,int,int) -> [str]
    Gibbs sampler
    """
    # get a list of randomly chosen t motifs
    dnas_kmers = [all_kmers(str, k) for str in dnas]
    motifs = [random.choice(k_mers) for k_mers in dnas_kmers]
    best_motifs = motifs # set above motifs to best motifs temporarily
    for j in range(n):
        index = random.choice(range(t))                                  # select an randomly generated index in range(0, t)
        motifs_subset = [motifs[i] for i in range(t) if i != index]      # remove an motif from motifs at randomly chosen index above
        profile = get_profile(motifs_subset)                             # get profile (:: [[probability]])
        motifs[index] = profile_randomly_gen_kmer(profile, dnas[index])  # reassign motif at the specific index chosen above
        if score(motifs) < score(best_motifs):                           # find best motifs by scoring
            best_motifs = motifs
    return best_motifs


def repeat(dnas, k, t, n, niter=20):
    """
    ([str],int,int,int,int) -> [str]
    """
    best_motifs = []
    best_score = 2**16
    for j in range(niter):
        motifs = gibbs_sampler(dnas, k, t, n)
        cur_score = score(motifs)
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
    return best_motifs


def main():
    f = open('/home/wsl/rosalind/data/ba02g.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    numbers = [int(frag.strip()) for frag in lines[0].split()]
    k = numbers[0]
    t = numbers[1]
    N = numbers[2]
    dnas = lines[1:]
    f.close()

    start_time = time.time()
    for motif in repeat(dnas, k, t, N):
        print(motif)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
