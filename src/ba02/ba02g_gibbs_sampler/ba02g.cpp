/*
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
- Modern way to filter STL container?
  https://stackoverflow.com/questions/21204676/modern-way-to-filter-stl-container
- Weighted random numbers
  https://stackoverflow.com/questions/1761626/weighted-random-numbers
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

int get_score(vector<string>);
unordered_map<char, vector<double>> get_profile_pseudocount(vector<string>, int);
double get_probability(string, unordered_map<char, vector<double>>);
string profile_randomly_generated_kmer(string, unordered_map<char, vector<double>>);
vector<string> gibbs_sampler(vector<string>, int, int, int);
vector<string> repeat_n_times(vector<string>, int, int, int, int);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02g.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::vector<std::string> splitted = split_str(lines[0], " ");
    int k = std::stoi(splitted[0]);
    int t = std::stoi(splitted[1]);
    int n = std::stoi(splitted[2]);
    std::vector<std::string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    clock_t tStart = clock();
    for (auto motif : repeat_n_times(dnas, k, t, n, 20)) {
        cout << motif << endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}

/************************************************
definitions
************************************************/

// 'Score(Motifs)' function
int get_score(vector<string> motifs) {
    const int t = motifs.size();
    const int k = motifs[0].length();
    int score = 0;
    for (int col = 0; col < k; col++) {
        unordered_map<char, double> temp = { {'A',0},{'C',0},{'G',0},{'T',0} };
        for (int row = 0; row < t; row++) {
            switch (motifs[row][col]) {
            case 'A':
                temp['A']++;
                break;
            case 'C':
                temp['C']++;
                break;
            case 'G':
                temp['G']++;
                break;
            case 'T':
                temp['T']++;
                break;
            }
        }
        int max_count = 0;
        for (auto const& p : temp) {
            if (p.second > max_count) {
                max_count = p.second;
            }
        }
        score += (t - max_count);
    }
    return score;
}

// 'Profile(Motifs)' function
unordered_map<char, vector<double>> get_profile_pseudocount(vector<string> motifs, int pseudocount) {
    unordered_map<char, vector<double>> profile = { {'A',{}},{'C',{}},{'G',{}},{'T',{}} };
    const int t = motifs.size();
    const int k = motifs[0].length();
    for (int col = 0; col < k; col++) {
        unordered_map<char, double> temp = { {'A',0.0},{'C',0.0},{'G',0.0},{'T',0.0} };
        for (int row = 0; row < t; row++) {
            switch (motifs[row][col]) {
            case 'A':
                temp['A']++;
                break;
            case 'C':
                temp['C']++;
                break;
            case 'G':
                temp['G']++;
                break;
            case 'T':
                temp['T']++;
                break;
            }
        }
        for (auto it = temp.begin(); it != temp.end(); it++) {
            profile[it->first].push_back((it->second + pseudocount) / (double)(t + t * pseudocount));
        }
    }
    return profile;
}

// 'Pr(pattern, profile)' function
double get_probability(string pattern, unordered_map<char, vector<double>> profile) {
    double prob = 1.0;
    for (size_t i = 0; i < pattern.length(); i++) {
        prob *= profile[pattern[i]][i];
    }
    return prob;
}

// profile randomly generated k-mer
// https://stackoverflow.com/questions/1761626/weighted-random-numbers
string profile_randomly_generated_kmer(string dna, unordered_map<char, vector<double>> profile) {
    int k = profile['A'].size();
    vector<double> weights;
    vector<string> kmers;
    for (size_t i = 0; i < dna.length() - k + 1; i++) {
        string kmer = dna.substr(i, k);
        kmers.push_back(kmer);
        weights.push_back(get_probability(kmer, profile));
    }
    std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
    std::mt19937 gen;
    gen.seed(time(0));
    return kmers[dist(gen)];
}

// BA02G: Gibbs sampler
vector<string> gibbs_sampler(vector<string> dnas, int k, int t, int n) {
    srand(time(NULL));
    vector<string> motifs;
    for (int i = 0; i < t; i++) {
        int idx = rand() % t;
        motifs.push_back(dnas[i].substr(idx, k));
    }
    vector<string> best_motifs = motifs;
    for (int i = 1; i < n; i++) {
        int idx = rand() % t;
        vector<string> new_motifs;
        int k = 0;
        while (k < t) {
            if (k == idx) {
                k++;
                continue;
            }
            new_motifs.push_back(motifs[k++]);
        }
        unordered_map<char, vector<double>> profile = get_profile_pseudocount(new_motifs, 1);
        motifs[idx] = profile_randomly_generated_kmer(dnas[idx], profile);
        if (get_score(motifs) < get_score(best_motifs)) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

// repeat n times
// not very fast
vector<string> repeat_n_times(vector<string> dnas, int k, int t, int n, int repeat) {
    vector<string> best_motifs;
    int best_score = std::numeric_limits<int>::max();
    for (int i = 0; i < repeat; i++) {
        vector<string> motifs = gibbs_sampler(dnas, k, t, n);
        int score = get_score(motifs);
        if (score < best_score) {
            best_score = score;
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

