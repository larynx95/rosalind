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
      -> GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> PREV: Randomized Motif Search (BA2F)
      -> HERE: Gibbs Sampling (BA2G)

Info:

    ttaccttAAC    tTACcttaac            ttaccttAAC    ttaccttAAC
    gATAtctgtc    gatATCtgtc            gATAtctgtc    gatatcTGTc
    ACGgcgttcg -> acggcgTTCg            ACGgcgttcg -> ACGgcgttcg
    ccctAAAgag    ccctaaAGAg            ccctAAAgag    ccctAAAgag
    cgtcAGAggt    CGTcagaggt            cgtcAGAggt    cgtcAGAggt

    RandomizedMotifSearch               GibbsSampler
    (may change all k-mers in one step) (change one-k-mer in one step)

  (1) first process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ttACCTtaac    taac     A: 2 1 1 1     A: 2/4 1/4 1/4 1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 1 1 1     C: 0   1/4 1/4 1/4
  Dna ccgGCGTtag  -> ---------- ->  /    -> G: 1 1 1 0  -> G: 1/4 1/4 1/4 0
      cactaACGAg     cactaACGAg    acta     T: 1 1 1 2     T: 1/4 1/4 1/4 2/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 2 2 2     A: 3/8 2/8 2/8 2/8
                      adding pseudocount -> C: 1 2 2 2  -> C: 1/8 2/8 2/8 2/8
                                            G: 2 2 2 1     G: 2/8 2/8 2/8 1/8
                                            T: 2 2 2 3     T: 2/8 2/8 2/8 3/8

    k-mers from deleted DNA string: ccgGCGTtag
    ccgG    cgGC    gGCG    GCGT    CGTt    GTta    Ttag
    0       0       0       1/128   0       1/256   0

    k-mers from deleted DNA string: ccgGCGTtag
    ccgG    cgGC    gGCG    *GCGT*  CGTt    GTta    Ttag
    4/8^4   8/8^4   8/8^4   24/8^4  12/8^4  16/8^4  8/8^4  total=80/8^4
    ─────   ─────   ─────   ──────  ──────  ──────  ─────
    80/8^4  80/8^4  80/8^4  80/8^4  80/8^4  80/8^4  80/8^4
=>  4/80    8/80    8/80    24/80   12/80   16/80   8/80   divided by total=80/8^4
=> RANDOM(4/80, 8/80, 8/80, 24/80, 12/80, 16/80, 8/80)     hypothetical seven-sided die

  (2) second process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ----------     /       A: 2 0 0 1     A: 2/4 0   0   1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 2 1 0     C: 0   2/4 1/4 0
  Dna ccgGCGTtag  -> ccgGCGTtag -> GCGT* -> G: 2 1 2 0  -> G: 2/4 1/4 2/4 0
      cactaACGAg     cactaACGAg    acta     T: 0 1 1 3     T: 0   1/4 1/4 3/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 1 1 2     A: 3/8 1/8 1/8 2/8
                      adding pseudocount -> C: 1 3 2 1  -> C: 1/8 3/8 2/8 1/8
                                            G: 3 2 3 1     G: 3/8 2/8 3/8 1/8
                                            T: 1 2 2 4     T: 1/8 2/8 2/8 4/8
    k-mers from deleted DNA string: ttACCTtaac
    ttAC    tACC    *ACCT*  CCTt    CTta    Ttaa    taac
    2/8^4   2/8^4   72/8^4  24/8^4  8/8^4   4/8^4   1/8^4

  (3) third process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ttACCTtaac    ACCT*    A: 2 0 0 1     A: 2/4 0   0   1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 2 1 0     C: 0   2/4 1/4 0
  Dna ccgGCGTtag  -> ccgGCGTtag -> GCGT  -> G: 2 1 2 0  -> G: 2/4 1/4 2/4 0
      cactaACGAg     ----------     /       T: 0 1 1 3     T: 0   1/4 1/4 3/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 1 1 2     A: 3/8 1/8 1/8 2/8
                      adding pseudocount -> C: 1 3 2 1  -> C: 1/8 3/8 2/8 1/8
                                            G: 3 2 3 1     G: 3/8 2/8 3/8 1/8
                                            T: 1 2 2 4     T: 1/8 2/8 2/8 4/8
    k-mers from deleted DNA string: cactaACGAg
    cact    acta    ctaA    taAC    aACG    *ACGA*  CGAg
    15/8^4  9/8^4   2/8^4   1/8^4   9/8^4   27/8^4  2/8^4


═════════════════════════════════════════════════

References:
- Cumulative distribution function (CDF)
  https://en.wikipedia.org/wiki/Cumulative_distribution_function
- A weighted version of random.choice
  https://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
- get random integer from probability distribution
  Generate random numbers with a given (numerical) distribution
  https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution
- Are tuples more efficient than lists in Python?
  https://stackoverflow.com/questions/68630/are-tuples-more-efficient-than-lists-in-python
"""

import random
import time


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
    motif = random.choices(kmers, weights=norm_probs, k=1)  # <-- Notice this line.
    return motif.pop()


def gibbs_sampler(dnas, k, t, n):
    """
    ([str],int,int,int) -> [str]
    implementation of 'GibbsSampler' function
    """
    motifs = random_motifs(dnas, k)
    best_motifs = motifs
    for j in range(n):
        index = random.choice(range(t))
        profile = get_profile([motifs[i] for i in range(t) if i != index])
        motifs[index] = profile_randomly_gen_kmer(profile, dnas[index])
        if score(motifs) < score(best_motifs):
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
    n = numbers[2]
    dnas = lines[1:]
    f.close()

    start_time = time.time()
    for motif in repeat(dnas, k, t, n):
        print(motif)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
