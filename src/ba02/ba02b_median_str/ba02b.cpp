/*
Rosalind: BA2A
Find a Median String

Given a k-mer Pattern and a longer string Text,
we use d(Pattern, Text) to denote
the minimum Hamming distance between Pattern and any k-mer in Text,

d(Pattern, Text) = min                 HammingDistance(Pattern, Pattern')
                  all k-mers Pattern' in Text

Given a k-mer Pattern and a set of strings
Dna = {Dna_1, ... , Dna_t},
we define d(Pattern, Dna) as the sum of distances
between Pattern and all strings in Dna,

                  t
d(Pattern, Dna) = Σ  d(Pattern, Dna_i)
                 i=0

Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna) over all
k-mers Pattern, the same task that the Equivalent Motif Finding Problem is
trying to achieve. We call such a k-mer a median string for Dna.

Median String Problem
Find a median string.

Given: An integer k and a collection of strings Dna.

Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern.
(If multiple answers exist, you may return any one.)

Sample Dataset
3
AAATTGACGCAT
GACGACCACGTT
CGTCAGCGCCTG
GCTGAGCACCGG
AGTACGGGACAG

Sample Output
GAC

═════════════════════════════════════════════════
References:
- Return an empty set with "return std::set<int>()" - why does it run?
  https://stackoverflow.com/questions/25856451/return-an-empty-set-with-return-stdsetint-why-does-it-run
  std::set<int>()

*/

#include <limits>
#include <cmath>
#include <unordered_set>
#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

int get_hamming_distance(std::string const& str1, std::string const& str2);
std::string number_to_pattern(int number, int num_k);
int distance_bw_pattern_strings(std::string const& pattern, std::vector<std::string> const& texts);
std::unordered_set<std::string> median_string(std::vector<std::string> const& dnas, int k);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02b.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    int k = std::stoi(lines[0]);
    std::vector<std::string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    // check time
    clock_t tStart = clock();
    for (auto kmer : median_string(dnas, k)) {
        std::cout << kmer << std::endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}

/************************************************
definitions
************************************************/

// BA01G: hamming distance
int get_hamming_distance(std::string const& str1, std::string const& str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
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

// BA02H: distance between pattern and strings
int distance_bw_pattern_strings(std::string const& pattern, std::vector<std::string> const& dnas) {
    int k = pattern.length();
    int distance = 0;
    for (auto dna : dnas) {
        int min_hdist = std::numeric_limits<int>::max(); // integer infinity
        for (int i = 0; i < dna.length() - k + 1;i++) {
            std::string kmer = dna.substr(i, k);
            int hdist = get_hamming_distance(kmer, pattern);
            if (hdist < min_hdist) {
                min_hdist = hdist;
            }
        }
        distance += min_hdist;
    }
    return distance;
}

// BA02H: median string
std::unordered_set<std::string> median_string(std::vector<std::string> const& dnas, int k) {
    int min_distance = std::numeric_limits<int>::max();
    std::unordered_set<std::string> median;
    for (int i = 0; i < std::pow(4, k);i++) {
        std::string kmer = number_to_pattern(i, k);
        int dist = distance_bw_pattern_strings(kmer, dnas);
        if (dist < min_distance) {
            min_distance = dist;
            median = std::unordered_set<std::string>();
            median.insert(kmer);
        }
        else if (dist == min_distance) {
            median.insert(kmer);
        }
    }
    return median;
}
