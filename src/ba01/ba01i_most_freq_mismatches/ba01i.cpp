/*
Rosalind: BA1I
Find the Most Frequent Words with Mismatches in a String

We defined a mismatch in "Compute the Hamming Distance Between Two Strings". We
now generalize "Find the Most Frequent Words in a String" to incorporate
mismatches as well.

Given strings Text and Pattern as well as an integer d,
we define Count_d(Text, Pattern)
as the total number of occurrences of Pattern in Text with at most d mismatches.
For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4
because AAAAA appears four times in this string with at most one mismatch:
AACAA, ATAAA, AAACA, and AAAGA.
Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in Text is simply a string Pattern
maximizing Count_d(Text, Pattern) among all k-mers. Note that Pattern does not
need to actually appear as a substring of Text; for example, AAAAA is the most
frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA
does not appear exactly in this string. Keep this in mind while solving the
following problem.

Frequent Words with Mismatches Problem
Find the most frequent k-mers with mismatches in a string.

Given: A string Text as well as integers k and d.

Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
GATG ATGC ATGT

═════════════════════════════════════════════════

Plan 1.
  * consider all possible 4^k k-mers, not just all k-mers (substrings) from a text

Pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  APPROXIMATEPATTERNCOUNT(Text, Pattern, d)        ║
  ║    count 0                                        ║
  ║    for i 0 to |Text| - |Pattern|                  ║
  ║      Pattern <- Text(i, |Pattern|)                ║
  ║      if HAMMINGDISTANCE(Pattern, Pattern’) >= d   ║
  ║        count count + 1                            ║
  ║    return count                                   ║
  ╚═══════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
  ║    Neighborhood <- the set consisting of the single string Pattern    ║
  ║    for i = 1 to |Pattern|                                             ║
  ║      symbol <- i-th nucleotide of Pattern                             ║
  ║      for each nucleotide x different from symbol                      ║
  ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
  ║        add Neighbor to Neighborhood                                   ║
  ║    return Neighborhood                                                ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NEIGHBORS(Pattern, d)                                   ║
  ║    if d = 0                                              ║
  ║      return {Pattern}                                    ║
  ║    if |Pattern| = 1                                      ║
  ║      return {A, C, G, T}                                 ║
  ║    Neighborhood <- an empty set                          ║
  ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
  ║    for each string Text from SuffixNeighbors             ║
  ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
  ║        for each nucleotide x                             ║
  ║          add x + Text to Neighborhood                    ║
  ║      else                                                ║
  ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
  ║    return Neighborhood                                   ║
  ╚══════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════╗
  ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
  ║    Neighborhood <- set consisting of single string Pattern   ║
  ║    for j = 1 to d                                            ║
  ║      for each string Pattern' in Neighborhood                ║
  ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
  ║        remove duplicates from Neighborhood                   ║
  ║    return Neighborhood                                       ║
  ╚══════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  FREQUENTWORDSWITHMISMATCHES(Text, k, d)                              ║
  ║    FrequentPatterns an empty set                                      ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLOSE(i) <- 0                                                    ║
  ║      FREQUENCYARRAY <- 0                                              ║
  ║    for i <- 0 to |Text| - k                                           ║
  ║      Neighborhood <- NEIGHBORS(Text(i, k), d)                         ║
  ║      for each Pattern from Neighborhood                               ║
  ║        index <- p2n(Pattern)                                          ║
  ║        CLOSE(index) <- 1                                              ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLOSE(i) = 1                                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        FREQUENCYARRAY(i) <- APPROXIMATEPATTERNCOUNT(Text, Pattern, d) ║
  ║    maxCount maximal value in FREQUENCYARRAY                           ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if FREQUENCYARRAY(i) = maxCount                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════════════════╗
  ║  FINDINGFREQUENTWORDSWITHMISMATCHESBYSORTING(Text, k, d)                 ║
  ║    FrequentPatterns <- an empty set                                      ║
  ║    Neighborhoods <- an empty list                                        ║
  ║    for i <- 0 to |Text| - k                                              ║
  ║      add NEIGHBORS(Text(i, k), d) to Neighborhoods                       ║
  ║    form an array NEIGHBORHOODARRAY holding all strings in Neighborhoods  ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      Pattern <- NEIGHBORHOODARRAY(i)                                     ║
  ║      INDEX(i) <-  p2n(Pattern)                                           ║
  ║      COUNT(i) <- 1                                                       ║
  ║    SORTEDINDEX SORT(INDEX)                                               ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      if SORTEDINDEX(i) = SORTEDINDEX(i + 1)                              ║
  ║        COUNT(i + 1) <- COUNT(i) + 1                                      ║
  ║    maxCount maximum value in array COUNT                                 ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║     if COUNT(i) = maxCount                                               ║
  ║       Pattern <- n2p(SORTEDINDEX(i), k)                                  ║
  ║       add Pattern to FrequentPatterns                                    ║
  ║    return FrequentPatterns                                               ║
  ╚══════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

Plan 1.
  (1) get all permutation with repetition for 4 nucleotides ('A', 'C', 'G', 'T')
      This is not efficient.
      TODO: Improve code with more efficient algorithm.
  (2) create an unordered_map
      key: each elem of (1)
      value: size of 'find_all_approx_occur(pattern, genome, d)'
  (3) get elem with max value in unordered_map

Plan 2.

Topic:
- permutation, combination
- iterator, generator
- recursion
- function: parameter type, decay, return type ...
*/

#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

int get_hamming_distance(std::string str1, std::string str2);
std::vector<int> find_all_approx_occur(std::string pattern, std::string genome,
    int num_d);
std::vector<std::string> perm_with_repetition(std::vector<std::string> vec,
    std::string nuc, int num);
std::vector<std::string> most_freq_words_with_mismatch(std::string genome,
    int num_k, int num_d);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01i.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string genome = lines[0];
    std::vector<std::string> numbers = split_str(lines[1], " ");
    int num_k = std::stoi(numbers[0]);
    int num_d = std::stoi(numbers[1]);

    clock_t tStart = clock();
    std::vector<std::string> answer =
        most_freq_words_with_mismatch(genome, num_k, num_d);
    print_vec(answer);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}

/************************************************
definitions
************************************************/

int get_hamming_distance(std::string str1, std::string str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

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

std::vector<std::string> most_freq_words_with_mismatch(std::string genome,
    int num_k, int num_d) {
    int max_freq = 0;
    std::unordered_map<std::string, int> dict;
    std::vector<std::string> kmers = perm_with_repetition({}, "ATGC", num_k);
    for (auto kmer : kmers) {
        int sz = find_all_approx_occur(kmer, genome, num_d).size();
        if (sz > max_freq) max_freq = sz;
        dict[kmer] = sz;
    }

    std::vector<std::string> answer;
    for (auto iter = dict.begin(); iter != dict.end(); iter++) {
        if (iter->second == max_freq) {
            answer.push_back(iter->first);
        }
    }
    return answer;
}