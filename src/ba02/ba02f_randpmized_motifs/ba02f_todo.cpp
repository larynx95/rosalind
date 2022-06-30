/*
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
* Random number generation
  - equivalent function to numpy.random.choice in C++
    https://stackoverflow.com/questions/42926209/equivalent-function-to-numpy-random-choice-in-c
  - Random number c++ in some range [duplicate]
    https://stackoverflow.com/questions/7560114/random-number-c-in-some-range
  - Generate random numbers uniformly over an entire range
    https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range
  - How to generate a random number in C++?
    https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c
*/

#include <random>
#include <unordered_set>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

int get_score(vector<string> motifs);
unordered_map<char, vector<double>> get_profile_pseudocount(vector<string> motifs);
double get_probability(string pattern, unordered_map<char, vector<double>> profile);
string profile_most_probable_kmer(string text, int k, unordered_map<char, vector<double>> profile);
vector<string> get_motifs_by_profile(vector<string> dnas, unordered_map<char, vector<double>> profile);
int get_random_integer(int n, int start_inclusive, int end_inclusive);
vector<string> get_random_motifs(vector<string> dnas, int k, int t);
vector<string> randomized_motif_search(vector<string> dnas, int k, int t);
vector<string> repeat_n_times(vector<string> dnas, int k, int t, int n);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02f.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::vector<std::string> splitted = split_str(lines[0], " ");
    int k = std::stoi(splitted[0]);
    int t = std::stoi(splitted[1]);
    std::vector<std::string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    clock_t tStart = clock();
    for (auto motif : repeat_n_times(dnas, k, t, 1000)) {
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
unordered_map<char, vector<double>> get_profile_pseudocount(vector<string> motifs, int pseudocount = 1) {
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

// BA02C: profile most probable k-mer
string profile_most_probable_kmer(string text, int k, unordered_map<char, vector<double>> profile) {
    string answer;
    double max_probability = -std::numeric_limits<double>::infinity();  // negative infinity
    for (size_t i = 0; i < text.length() - k + 1; i++) {
        string kmer = text.substr(i, k);
        double new_prob = get_probability(kmer, profile);
        if (new_prob > max_probability) {
            max_probability = new_prob;
            answer = kmer;
        }
    }
    return answer;
}

// 'Motifs(Profile, Dnas)', get motifs from profile
vector<string> get_motifs_by_profile(vector<string> dnas, unordered_map<char, vector<double>> profile) {
    vector<string> motifs;
    int t = dnas.size();
    int k = profile['A'].size();
    for (int i = 0; i < t; i++) {
        motifs.push_back(profile_most_probable_kmer(dnas[i], k, profile));
    }
    return motifs;
}

// get a random integer
int get_random_integer(int n, int start_inclusive, int end_inclusive) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(start_inclusive, end_inclusive);
    return distr(generator);
}

// get n random integers
// https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
vector<int> get_n_random_integers(int n, int start_inclusive, int end_inclusive) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(start_inclusive, end_inclusive);
    vector<int> ints;
    for (auto i = 0; i < n; i++) {
        ints.push_back(distr(generator));
    }
    return ints;
}

// get random motifs
vector<string> get_random_motifs(vector<string> dnas, int k, int t) {
    int n = dnas[0].length();
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(0, n - k);
    vector<string> motifs;
    for (int i = 0; i < t; i++) {
        int idx = distr(generator);
        motifs.push_back(dnas[i].substr(idx, k));
    }
    return motifs;
}

// BA02F: randomized motif search
// TODO: fix this function! no error, incorrect result
vector<string> randomized_motif_search(vector<string> dnas, int k, int t) {
    vector<string> motifs = get_random_motifs(dnas, k, t);
    vector<string> best_motifs = motifs;
    for (; true; ) {
        unordered_map<char, vector<double>> profile = get_profile_pseudocount(motifs, 1);
        vector<string> motifs = get_motifs_by_profile(dnas, profile);
        if (get_score(motifs) < get_score(best_motifs)) {
            best_motifs = motifs;
        }
        else {
            return best_motifs;
        }
    }
}

// repeat n times
vector<string> repeat_n_times(vector<string> dnas, int k, int t, int n) {
    int min_score = std::numeric_limits<int>::max();
    vector<string> best_motifs;
    for (int i = 0; i < n; i++) {
        vector<string> motifs = randomized_motif_search(dnas, k, t);
        int score = get_score(motifs);
        if (score < min_score) {
            min_score = score;
            best_motifs = motifs;
        }
    }
    return best_motifs;
}