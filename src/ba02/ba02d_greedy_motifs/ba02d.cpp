/*
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

Info.
                          T  C  G  G  G  G  g  T  T  T  t  t
                          c  C  G  G  t  G  A  c  T  T  a  C
                          a  C  G  G  G  G  A  T  T  T  t  C
                          T  t  G  G  G  G  A  c  T  T  t  t
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C
                          T  t  G  G  G  G  A  c  T  T  C  C
                          T  C  G  G  G  G  A  T  T  c  a  t
                          T  C  G  G  G  G  A  T  T  c  C  t
                          T  a  G  G  G  G  A  a  c  T  a  C
                          T  C  G  G  G  t  A  T  a  a  C  C

    SCORE(Motifs)         3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30

    COUNT(Motifs)     A:  2  2  0  0  0  0  9  1  1  1  3  0
                      C:  1  6  0  0  0  0  0  4  1  2  4  6
                      G:  0  0 10 10  9  9  1  0  0  0  0  0
                      T:  7  2  0  0  1  1  0  5  8  7  3  4

    PROFILE(Motifs)   A: .2 .2  0  0  0  0 .9 .1 .1 .1 .3  0
                      C: .1 .6  0  0  0  0  0 .4 .1 .2 .4 .6
                      G:  0  0  1  1 .9 .9 .1  0  0  0  0  0
                      T: .7 .2  0  0 .1 .1  0 .5 .8 .7 .3 .4
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

                      A: .2 .2 .0 .0 .0 .0 .9 .1 .1 .1 .3 .0
                      C: .1 .6 .0 .0 .0 .0 .0 .4 .1 .2 .4 .6
                      G: .0 .0  1  1 .9 .9 .1 .0 .0 .0 .0 .0
                      T: .7 .2 .0 .0 .1 .1 .0 .5 .8 .7 .3 .4
    Pr(ACGGGGATTACC|Profile) = .2*.6*1*1*.9*.9*.9*.5*.8*.1*.4*.6 = 0.000839808
    (profile most probable kmer is not the same as concensus!!)

Plan 1.
-
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

- TODO: GreedyMotifSearch fx in BA02D is not always correct.
  Because 'Profile-most probable k-mer' fx can have more than one result.
  To pass this BA02D, take the first k-mer in 'Profile-most probable k-mer',
  not select one from the result HashSet.

═════════════════════════════════════════════════

References:
* Greedy Algorithm
  https://www.programiz.com/dsa/greedy-algorithm#:~:text=A%20greedy%20algorithm%20is%20an,if%20the%20choice%20is%20wrong.
* modifying vector value in HashMap
  - How can I update a value in a mutable HashMap?
    https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap
  - HashMap or_insert()
    https://doc.rust-lang.org/std/collections/hash_map/enum.Entry.html#method.or_insert
  - How to build a HashMap of Vectors in Rust? (uncomfortable)
    https://stackoverflow.com/questions/26169216/how-to-build-a-hashmap-of-vectors-in-rust
- What is the difference between iter and into_iter?
  https://stackoverflow.com/questions/34733811/what-is-the-difference-between-iter-and-into-iter
- Using map with Vectors
  https://stackoverflow.com/questions/30026893/using-map-with-vectors
- If only one element in a hashset, how can I get it out?
  https://stackoverflow.com/questions/23595749/if-only-one-element-in-a-hashset-how-can-i-get-it-out
  set.iter().next().unwrap();
*/

#include <unordered_set>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

unordered_map<char, vector<int>> get_count(vector<string>);
int get_score(vector<string>);
unordered_map<char, vector<double>> get_profile(vector<string>);
double get_probability(string, unordered_map<char, vector<double>>);
string profile_most_probable_kmer(string, int, unordered_map<char, vector<double>>);
vector<string> greedy_motif_search(vector<string>, int, int);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02d.txt" : argv[1];
    vector<string> lines = read_lines(filename);
    vector<string> splitted = split_str(lines[0], " ");
    int k = stoi(splitted[0]);
    int t = stoi(splitted[1]);
    vector<string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    // check time
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
unordered_map<char, vector<double>> get_profile(vector<string> motifs) {
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
            profile[it->first].push_back(it->second / t);
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
            unordered_map<char, vector<double>> profile = get_profile(motifs);
            motifs.push_back(profile_most_probable_kmer(dnas[j], k, profile));
        }
        if (get_score(motifs) < get_score(best_motifs)) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}
