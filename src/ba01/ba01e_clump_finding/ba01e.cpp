/*
Rosalind: BA1E
Find Patterns Forming Clumps in a String

TODO: Find more efficient algorithms!

Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger)
string Genome if there is an interval of Genome of length L in which Pattern
appears at least t times. For example, TGCA forms a (25,3)-clump in the
following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem
Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset
CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4

Sample Output
CGACA GAAGA AATGT

═════════════════════════════════════════════════

Pseudocode:
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  CLUMPFINDING(Genome, k, t, L)                                        ║
  ║    FrequentPatterns <- an empty set                                   ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLUMP(i) <- 0                                                    ║
  ║    for i <- 0 to |Genome| - L                                         ║
  ║      Text <- the string of length L starting at position i in Genome  ║
  ║      FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)                  ║
  ║      for index <- 0 to 4k - 1                                         ║
  ║        if FREQUENCYARRAY(index) >= t                                  ║
  ║          CLUMP(index) <- 1                                            ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLUMP(i) = 1                                                  ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                               ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════╗
  ║  CLUMPFINDINGBetter(Genome, k, t, L)                      ║
  ║    FrequentPatterns <- an empty set                       ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      CLUMP(i) <- 0                                        ║
  ║    Text <- Genome(0, L)                                   ║
  ║    FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)        ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if FREQUENCYARRAY(i) - t                             ║
  ║        CLUMP(i) <- 1                                      ║
  ║    for i <- 1 to |Genome| - L                             ║
  ║      FirstPattern <- Genome(i - 1, k)                     ║
  ║      index <- PATTERNTONUMBER(FirstPattern)               ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) - 1   ║
  ║      LastPattern <- Genome(i + L - k, k)                  ║
  ║      index <- PATTERNTONUMBER(LastPattern)                ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) + 1   ║
  ║      if FREQUENCYARRAY(index) >= t                        ║
  ║        CLUMP(index) <- 1                                  ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if CLUMP(i) = 1                                      ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                   ║
  ║        add Pattern to the set FrequentPatterns            ║
  ║    return FrequentPatterns                                ║
  ╚═══════════════════════════════════════════════════════════╝

1. find the generalized rule of frequency array

  "AAGCAAAGGTGGGC"  len=14
   vv-----------                                    first ... last
  "AAGCAAAGGTGGG"   len=13, starting from index 0   AA    ...  -
   "AGCAAAGGTGGGC"  len=13, starting from index 1   -     ...  GC
    -----------^^
    common part!

  (1) compute_freq("AAGCAAAGGTGGG", 2)
                    ^^  <--- minus 1
  (2) compute_freq( "AGCAAAGGTGGGC", 2)
                                ^^  <--- plus 1

      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
      3  0  2  0  1  0  0  0  0  1  3  1  0  0  1  0   --- (1)
     -2  0  2  0  1  0  0  0  0  2+ 3  1  0  0  1  0   --- (2)

2. If we know a frequency array of the first clump,
   we can get the frequency array of the whole genome.

    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT AAGCAAAGGTGGGC   fst  lst
    3  0  2  0  1  0  0  0  0  1  1  1  0  0  1  0  AAGCAAAGGTG
   -2  0  2  0  1  0  0  0  0  1  2+ 1  0  0  1  0   AGCAAAGGTGG     AA   GG
    2  0 -1  0  1  0  0  0  0  1  3+ 1  0  0  1  0    GCAAAGGTGGG    AG   GG
    2  0  1  0  1  0  0  0  0 -1+ 3  1  0  0  1  0     CAAAGGTGGGC   GC   GC
    check the frequency before subtraction!

═════════════════════════════════════════════════

Plan 1:
  - given strings:
    strings   n
    ------------------------
    DNA       just one string
    clumps    (len(DNA) - L) for DNA
    k-mers    (L-k) for each L-clump
  - For each L-clump, find the frequencies for k-mers in the L-clump.
  - Repeat above process for all L-clumps.

Plan 2:
  - using 'frequency array'?
  - TODO: Find efficient algorithm!

Topics:
  1. C++, language specific topics
  - split string
  - string to int conversion
  - check thether a vector contains an element
  - array push_back
  2. algorithmic thinking - about better algorithm
  - select faster compiling language: C, C++
  - write down plans before typing keyboard
*/

#include <time.h>
#include <cmath>
#include <set>
#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

template <typename T>
bool is_vec_contains(std::vector<T> vec, T elem);

std::vector<int> find_all_ocurrences(std::string pattern, std::string genome);
std::vector<std::string> clump_finding_perm(const std::string genome, int num_k,
    int num_L, int num_t);

long long int pattern_to_number(std::string);
std::string number_to_pattern(int, int);
std::vector<int> compute_freq(std::string, int);
std::set<std::string> clump_finding(const std::string genome, int num_k,
    int num_L, int num_t);
std::set<std::string> clump_finding_better(const std::string genome, int num_k,
    int num_L, int num_t);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01e.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string genome = lines[0];
    std::string numbers = lines[1];
    std::vector<std::string> splitted = split_str(numbers, " ");  // split_str fx.
    int num_k = std::stoi(splitted[0]);
    int num_L = std::stoi(splitted[1]);
    int num_t = std::stoi(splitted[2]);
    int num_genome = genome.length();

    /*
    clock_t tStart1 = clock();
    std::vector<std::string> distinct_kmers_perm =
        clump_finding_perm(genome, num_k, num_L, num_t);
    print_vec(distinct_kmers_perm);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart1) / CLOCKS_PER_SEC);
    // Time taken: 77.73s

    clock_t tStart2 = clock();
    std::set<std::string> distinct_kmers =
        clump_finding(genome, num_k, num_L, num_t);
    print_set(distinct_kmers);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart2) / CLOCKS_PER_SEC);
    // Time taken: 87.73s
    */

    clock_t tStart3 = clock();
    std::set<std::string> distinct_kmers_better =
        clump_finding_better(genome, num_k, num_L, num_t);
    print_set(distinct_kmers_better);
    printf("Time taken: %.6fs\n", (double)(clock() - tStart3) / CLOCKS_PER_SEC);
    // Time taken: 0.06s

    return 0;
}

/************************************************
definitions
************************************************/

template <typename T>
bool is_vec_contains(std::vector<T> vec, T elem) {
    for (auto e : vec) {
        if (e == elem) return true;
    }
    return false;
}

// BA01D: find all occurrences
std::vector<int> find_all_ocurrences(std::string pattern, std::string genome) {
    std::vector<int> indices;
    int len_genome = genome.length();
    int len_pattern = pattern.length();

    for (int i = 0; i < len_genome - len_pattern + 1; i++) {
        if (genome.substr(i, len_pattern) == pattern) {
            indices.push_back(i);
        }
    }
    return indices;
}

// BA01E: clump finding
std::vector<std::string> clump_finding_perm(const std::string genome, int num_k,
    int num_L, int num_t) {
    std::vector<std::string> distinct_kmers;

    for (int i = 0; i < genome.size() - num_L + 1; i++) {
        std::string clump = genome.substr(i, num_L);
        for (int j = 0; j < clump.length() - num_k + 1; j++) {
            std::string kmer = clump.substr(j, num_k);
            int freq = find_all_ocurrences(kmer, clump).size();
            if (freq >= num_t && !is_vec_contains(distinct_kmers, kmer)) {
                distinct_kmers.push_back(kmer);
            }
        }
    }
    return distinct_kmers;
}

// BA01L: pattern to number
long long int pattern_to_number(std::string pattern) {
    long long int number = 0;
    int len = pattern.size();
    for (int i = 0; i < len; i++) {
        int pw = len - i - 1;
        int nuc;
        switch (pattern[i]) {
        case 'A':
            nuc = 0;
            break;
        case 'C':
            nuc = 1;
            break;
        case 'G':
            nuc = 2;
            break;
        case 'T':
            nuc = 3;
            break;
        }
        long double pw_val = std::pow(4.0, pw * 1.0);
        number += (pw_val * nuc);
    }
    return number;
}

// BA01M: number to pattern
std::string number_to_pattern(int number, int num_k) {
    std::string pattern;
    for (int i = 0; i < num_k; i++) {
        int temp = number / std::pow(4, num_k - i - 1);
        switch (temp) {
        case 0:
            pattern += 'A';
            break;
        case 1:
            pattern += 'C';
            break;
        case 2:
            pattern += 'G';
            break;
        case 3:
            pattern += 'T';
            break;
        }
        number -= temp * std::pow(4, num_k - i - 1);
    }
    return pattern;
}

// BA01K: compute frequency
std::vector<int> compute_freq(std::string pattern, int num_k) {
    std::vector<int> freq_vec(pow(4, num_k), 0);
    for (int i = 0; i < pattern.size() - num_k + 1; i++) {
        std::string kmer = pattern.substr(i, num_k);
        int temp = pattern_to_number(kmer);
        freq_vec[temp] += 1;
    }
    return freq_vec;
}

// BA01E: clump finding
std::set<std::string> clump_finding(const std::string genome, int k, int L,
    int t) {
    std::set<std::string> freq_patterns;
    std::vector<int> clump_freq_vec(pow(4, k), 0);
    for (int i = 0; i < genome.size() - L + 1; i++) {
        std::string clump = genome.substr(i, L);
        std::vector<int> freq_vec = compute_freq(clump, k);
        for (int i = 0; i < freq_vec.size(); i++) {
            if (freq_vec[i] >= t) clump_freq_vec[i] = 1;
        }
    }
    for (int i = 0; i < clump_freq_vec.size(); i++) {
        if (clump_freq_vec[i] == 1) {
            std::string pattern = number_to_pattern(i, k);
            freq_patterns.insert(pattern);
        }
    }
    return freq_patterns;
}

// BA01E: clump finding, better
std::set<std::string> clump_finding_better(const std::string genome, int k,
    int L, int t) {
    std::set<std::string> patterns;
    std::vector<int> clump_freq_vec(pow(4, k), 0);  // size 4^k vector, all 0
    std::string fst_clump = genome.substr(0, L);
    std::vector<int> freq_vec = compute_freq(fst_clump, k);
    for (int i = 1; i < genome.size() - L + 1; i++) {
        std::string fst_kmer = genome.substr(i - 1, k);
        int idx = pattern_to_number(fst_kmer);
        freq_vec[idx]--;
        std::string lst_kmer = genome.substr(i + L - k, k);
        idx = pattern_to_number(lst_kmer);
        freq_vec[idx]++;
        if (freq_vec[idx] >= t) {
            clump_freq_vec[idx] = 1;
        }
    }
    for (int i = 0; i < (int)pow(4, k); i++) {
        if (clump_freq_vec[i] == 1) {
            std::string pattern = number_to_pattern(i, k);
            patterns.insert(pattern);
        }
    }
    return patterns;
}
