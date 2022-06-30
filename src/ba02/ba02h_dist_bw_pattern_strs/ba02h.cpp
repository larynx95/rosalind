/*
Rosalind: BA2H
Implement DistanceBetweenPatternAndStrings

The first potential issue with implementing MedianString from "Find a Median String" is
writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai),
the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}.
This task is achieved by the following pseudocode.

  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

Sample Dataset
AAA
TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

Sample Output
5

═════════════════════════════════════════════════
Pseudocode:
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ MEDIANSTRING(Dna, k)                                                   ║
  ║     distance 1                                                         ║
  ║     for i 0 to 4^k - 1                                                 ║
  ║         pattern <- NumberToPattern(i, k)                               ║
  ║         if distance > DistanceBetweenPaternAndStrings(Pattern, Dna)    ║
  ║             distance <- DistanceBetweenPatternAndStrings(Patern, Dna)  ║
  ║             Median <- Pattern                                          ║
  ║     return Median                                                      ║
  ╚════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════
References:
- Passing std::string by Value or Reference [duplicate]
  https://stackoverflow.com/questions/10789740/passing-stdstring-by-value-or-reference
  1) Using the string as an id (will not be modified).
     Passing it in by const reference is probably the best idea here: (std::string const&)
  2) Modifying the string but not wanting the caller to see that change.
     Passing it in by value is preferable: (std::string)
  3) Modifying the string but wanting the caller to see that change.
     Passing it in by reference is preferable: (std::string &)
  4) Sending the string into the function and the caller of the function will never use the string again.
     Using move semantics might be an option (std::string &&)
- std::pow, std::powf, std::powl
  https://en.cppreference.com/w/cpp/numeric/math/pow
*/

#include <limits>
#include <unordered_set>
#include <cmath>  // pow
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
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02h.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    std::vector<std::string> dnas = split_str(lines[1], " ");

    // check time, distance between pattern and strings
    clock_t tStart = clock();
    std::cout << distance_bw_pattern_strings(pattern, dnas) << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // check time, median string
    tStart = clock();
    for (auto kmer : median_string(dnas, 3)) {
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
    for (size_t i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

// BA01M: number to pattern
std::string number_to_pattern(int number, int num_k) {
    std::string pattern;
    for (size_t i = 0; i < (size_t)num_k; i++) {
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
        for (size_t i = 0; i < dna.length() - k + 1;i++) {
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
