/*
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

═════════════════════════════════════════════════
References:
-
*/

#include <unordered_set>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

unordered_map<char, vector<int>> get_count(vector<string>);
int get_score(vector<string>);
unordered_map<char, vector<double>> get_profile_pseudocount(vector<string> motifs);  // <-- modified just this function
double get_probability(string, unordered_map<char, vector<double>>);
string profile_most_probable_kmer(string, int, unordered_map<char, vector<double>>);
vector<string> greedy_motif_search(vector<string>, int, int);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02e.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::vector<std::string> splitted = split_str(lines[0], " ");
    int k = std::stoi(splitted[0]);
    int t = std::stoi(splitted[1]);
    std::vector<std::string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    clock_t tStart = clock();
    for (auto motif : greedy_motif_search(dnas, k, t)) {
        cout << motif << endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}

/************************************************
definitions
************************************************/

// 'Count(Motifs)' function
unordered_map<char, vector<int>> get_count(vector<string> motifs) {
    unordered_map<char, vector<int>> count = { {'A',{}},{'C',{}},{'G',{}},{'T',{}} };
    const int t = motifs.size();
    const int k = motifs[0].length();
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
        for (auto nuc : "ACGT") {
            count[nuc].push_back(temp[nuc]);
        }
    }
    return count;
}

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

// BA02D: greedy motif search
vector<string> greedy_motif_search(vector<string> dnas, int k, int t) {
    vector<string> best_motifs;
    for (auto dna : dnas) {
        best_motifs.push_back(dna.substr(0, k));
    }
    for (size_t i = 0; i < dnas[0].length() - k + 1; i++) {
        vector<string> motifs = { dnas[0].substr(i, k) };
        for (int j = 1; j < t; j++) {
            unordered_map<char, vector<double>> profile = get_profile_pseudocount(motifs, 1);
            motifs.push_back(profile_most_probable_kmer(dnas[j], k, profile));
        }
        if (get_score(motifs) < get_score(best_motifs)) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}
