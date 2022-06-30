/*
Rosalind: BA1K
Generate the Frequency Array of a String

Given an integer k, we define the frequency array of a string Text as an array
of length 4k, where the i-th element of the array holds the number of times that
the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

kmer      AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
index      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
frequency  3   0   2   0   1   0   0   0   0   1   3   1   0   0   1   0

Computing a Frequency Array
Generate the frequency array of a DNA string.

Given: A DNA string Text and an integer k.

Return: The frequency array of k-mers in Text.

Sample Dataset
ACGCGGCTCTGAAA
2

Sample Output
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

═════════════════════════════════════════════════
Pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  COMPUTINGFREQUENCIES(Text, k)                    ║
  ║    for i <- 0 to 4^k- 1                           ║
  ║      FREQUENCYARRAY(i) <- 0                       ║
  ║    for i <- 0 to |Text| - k                       ║
  ║      Pattern <- Text(i, k)                        ║
  ║      j <- PATTERNTONUMBER(Pattern)                ║
  ║      FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1   ║
  ║    return FREQUENCYARRAY                          ║
  ╚═══════════════════════════════════════════════════╝

═════════════════════════════════════════════════
Plan 1.
  (1) get all k-mer patterns: permutation with repetition
  (2) sort k-mer patterns above
  (3) define a function for calculating frequency
      - 'std::count'
  (4) get frequency array

Plan 2.
  (1) define helper functions:
      - number to pattern
      - pattern to number
  (2) create an empty array with its size 4^n
  (3) figure out the relationship between the index and the pattern

Topic:
-
*/

#include <cmath>
#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
function prototypes
************************************************/

template <typename T>
void print_vec(const T& t);

int get_num_lines(std::string filename);
std::vector<std::string> get_lines(std::string filename);
std::vector<std::string> split_str(std::string str, std::string delimiter);
std::vector<std::string> perm_with_repetition(std::vector<std::string> vec,
    std::string nuc, int num);
int get_freq(const std::string&, const std::string&);
std::vector<int> compute_freq_brute(std::string, int);

long long int pattern_to_number(std::string);
std::string number_to_pattern(int, int);
std::vector<int> compute_freq(std::string, int);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01k.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    int num_k = std::stoi(lines[1]);

    std::vector<int> freq_vec = compute_freq(pattern, num_k);
    print_vec(freq_vec);

    return 0;
}

/************************************************
function definitions
************************************************/

std::vector<std::string> perm_with_repetition(std::vector<std::string> vec,
    std::string nuc, int num) {
    if (num == 0) {
        return vec;
    }
    else if (vec.size() == 0) {
        return perm_with_repetition({ "A", "C", "G", "T" }, nuc, num - 1);
    }
    else {
        std::vector<std::string> vec_new;
        for (auto elem : vec) {
            for (auto& chr : nuc) {
                std::string temp = elem + chr;
                vec_new.push_back(temp);
            }
        }
        vec = std::move(vec_new);
        return perm_with_repetition(vec, nuc, num - 1);
    }
}

int get_freq(const std::string& genome, const std::string& pattern) {
    int cnt = 0;
    for (int i = 0; i < genome.size() - pattern.size() + 1; i++) {
        if (pattern == genome.substr(i, pattern.size())) cnt++;
    }
    return cnt;
}

// BA01K: compute frequency, brute
std::vector<int> compute_freq_brute(std::string genome, int num_k) {
    std::vector<std::string> perm = perm_with_repetition({}, "ACGT", num_k);
    std::vector<int> freq_vec;  // TODO: problem here
    for (auto kmer : perm) {
        freq_vec.push_back(get_freq(genome, kmer));
    }
    return freq_vec;
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