/*
Rosalind: BA1B
Find the Most Frequent Words in a String

We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text,
Pattern) among all k-mers. For example, "ACTAT" is a most frequent 5-mer in
"ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most frequent 3-mer of
"CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4

Sample Output
CATG GCAT

═════════════════════════════════════════════════

Pseudocode:
╔═════════════════════════════════════════════════╗
║ COMPUTINGFREQUENCIES(Text, k)                   ║
║   for i <- 0 to 4^k- 1                          ║
║     FREQUENCYARRAY(i) <- 0                      ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- Text(i, k)                       ║
║     j <- PATTERNTONUMBER(Pattern)               ║
║     FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1  ║
║   return FREQUENCYARRAY                         ║
╚═════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════╗
║ FREQUENTWORDS(Text, k)                          ║
║   FrequentPatterns an empty set                 ║
║     for i <- 0 to |Text| - k                    ║
║     Pattern <- the k-mer Text(i, k)             ║
║     COUNT(i) <- PATTERNCOUNT(Text, Pattern)     ║
║   maxCount maximum value in array COUNT         ║
║   for i <- 0 to |Text| - k                      ║
║     if COUNT(i) = maxCount                      ║
║       add Text(i, k) to FrequentPatterns        ║
║   remove duplicates from FrequentPatterns       ║
║   return FrequentPatterns                       ║
╚═════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════╗
║ FASTERFREQUENTWORDS(Text , k)                     ║
║   FrequentPatterns <- an empty set                ║
║   FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k) ║
║   maxCount <- maximal value in FREQUENCYARRAY     ║
║   for i <- 0 to 4k - 1                            ║
║     if FREQUENCYARRAY(i) = maxCount               ║
║       Pattern <- NUMBERTOPATTERN(i, k)            ║
║       add Pattern to the set FrequentPatterns     ║
║   return FrequentPatterns                         ║
╚═══════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════╗
║ FINDINGFREQUENTWORDSBYSORTING(Text , k)              ║
║   FrequentPatterns <- an empty set                   ║
║   for i <- 0 to |Text| - k                           ║
║     Pattern <- Text(i, k)                            ║
║     INDEX(i) <- PATTERNTONUMBER(Pattern)             ║
║     COUNT(i) <- 1                                    ║
║   SORTEDINDEX <- SORT(INDEX)                         ║
║   for i <- 1 to |Text| - k                           ║
║     if SORTEDINDEX(i) = SORTEDINDEX(i - 1)           ║
║       COUNT(i) = COUNT(i - 1) + 1                    ║
║   maxCount <- maximum value in the array COUNT       ║
║   for i <- 0 to |Text| - k                           ║
║     if COUNT(i) = maxCount                           ║
║       Pattern <- NUMBERTOPATTERN(SORTEDINDEX(i), k)  ║
║       add Pattern to the set FrequentPatterns        ║
║   return FrequentPatterns                            ║
╚══════════════════════════════════════════════════════╝

Plan 1.
  - Extract k-mer in order from the beginning of the given pattern
    and record the frequency at the same time.
  - using two arrays or map (unordered_map)

═════════════════════════════════════════════════
References:
- How to find if a given key exists in a C++ std::map
  https://stackoverflow.com/questions/1939953/how-to-find-if-a-given-key-exists-in-a-c-stdmap
- std::unordered_map
  https://en.cppreference.com/w/cpp/container/unordered_map
- C++ Loop through Map
  https://stackoverflow.com/questions/26281979/c-loop-through-map
*/

#include <map>
#include <set>
#include "/home/wsl/rosalind/include/utils.hpp"

// BA01B: frequent words
std::set<std::string> freq_words(std::string text, int k) {
    std::unordered_map<std::string, int> map;
    std::set<std::string> words;
    // construct map, find the max frequency
    int max_freq = 0;
    for (int i = 0; i < text.length() - k + 1; i++) {
        std::string pattern = text.substr(i, k);
        if (map.find(pattern) == map.end()) {  // if absent
            map[pattern] = 1;
        }
        else {
            map[pattern] += 1;
            if (map[pattern] > max_freq) {
                max_freq = map[pattern];
            }
        }
    }
    // find words with max frequency
    for (auto pair : map) {
        if (pair.second == max_freq) {
            words.insert(pair.first);
        }
    }
    return words;
}

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01b.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string text = lines[0];
    int k = std::stoi(lines[1]);

    // check time
    clock_t tStart = clock();
    for (auto e : freq_words(text, k)) {
        std::cout << e << " ";
    }
    printf("\n");
    printf("Time taken: %.5fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}
