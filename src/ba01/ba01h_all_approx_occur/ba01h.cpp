/*
Rosalind: BA1H
Find All Approximate Occurrences of a Pattern in a String

We say that a k-mer Pattern appears as a substring of Text with at most d
mismatches if there is some k-mer substring Pattern' of Text having d or fewer
mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our
observation that a DnaA box may appear with slight variations leads to the
following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return: All starting positions where Pattern appears as a substring of Text with
at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3

Sample Output
6 7 26 27 78

═════════════════════════════════════════════════

Pseudocode:

APPROXIMATEPATTERNCOUNT(Text, Pattern, d)
  count <- 0
  for i <- 0 to |Text| � |Pattern|
    Pattern' <- Text(i, |Pattern|)
    if HAMMINGDISTANCE(Pattern, Pattern’) <= d
      count <- count + 1
  return count

*/

#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

int get_hamming_distance(std::string str1, std::string str2);
std::vector<int> find_all_approx_occur(std::string pattern, std::string genome,
    int num_d);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01h.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    std::string genome = lines[1];
    int num_d = std::stoi(lines[2]);
    print_vec(find_all_approx_occur(pattern, genome, num_d));
    return 0;
}

/************************************************
definitions
************************************************/

// BA01G: hamming distance
int get_hamming_distance(std::string str1, std::string str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

// BA01H: all approximate occurrences
std::vector<int> find_all_approx_occur(std::string pattern, std::string genome,
    int num_d) {
    std::vector<int> indices;
    for (int i = 0; i < genome.size() - pattern.size() + 1; i++) {
        std::string frag = genome.substr(i, pattern.size());
        if (get_hamming_distance(pattern, frag) <= num_d) {
            indices.push_back(i);
        }
    }
    return indices;
}