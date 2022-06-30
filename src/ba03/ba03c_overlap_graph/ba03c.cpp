/*
Rosalind: BA3C
Construct the Overlap Graph of a Collection of k-mers

In this chapter, we use the terms prefix and suffix
to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.

Given an arbitrary collection of k-mers Patterns,
we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern'
by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern').
The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

Overlap Graph Problem
Construct the overlap graph of a collection of k-mers.

Given:
A collection Patterns of k-mers.

Return:
The overlap graph Overlap(Patterns), in the form of an adjacency list.

Sample Dataset
ATGCG
GCATG
CATGC
AGGCA
GGCAT

Sample Output (1:1)
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * HERE: Construct OverapGraph (BA3C)
      -> NEXT: Reconstruct a String from its k-mer Composition (BA3H)

Plan 1.
  * sample dataset
    kmer     prefix suffix
    ----------------------
    ATGCG -> ATGC   TGCG
    GCATG -> GCAT   CATG
    CATGC -> CATG   ATGC
    AGGCA -> AGGC   GGCA
    GGCAT -> GGCA   GCAT

═════════════════════════════════════════════════

References:
- Determine if map contains a value for a key?
  https://stackoverflow.com/questions/3136520/determine-if-map-contains-a-value-for-a-key
- td::unordered_map::at
  https://cplusplus.com/reference/unordered_map/unordered_map/at/
- Iterate through unordered map C++
  https://stackoverflow.com/questions/22880431/iterate-through-unordered-map-c
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

// BA03C: overlap graph
unordered_map<string, vector<string>> overlap_graph(vector<string> patterns) {
    int k = patterns[0].length();
    unordered_map<string, vector<string>> dic;
    for (auto source : patterns) {
        for (auto target : patterns) {
            if (source.substr(1) == target.substr(0, k - 1)) {
                if (!dic.count(source)) { dic.insert({ source, {target} }); }
                else dic.at(source).push_back(target);
            }
        }
    }
    return dic;
}

// main
int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03c.txt");

    clock_t tStart = clock();
    for (const auto& [key, value] : overlap_graph(lines)) {
        cout << key << " -> ";
        for (auto val : value) {
            cout << val << " ";
        }
        cout << endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}