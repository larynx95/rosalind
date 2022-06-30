/*
Rosalind: BA3E (difficulty: 1/5)
Construct the "De Bruijn Graph" of a Collection of k-mers

Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times),
we define "CompositionGraph(Patterns)" as a graph with |Patterns| isolated edges.
Every edge is labeled by a k-mer from Patterns,
and the starting and ending nodes of an edge are labeled
by the prefix and suffix of the k-mer labeling that edge.
We then define the "de Bruijn graph" of Patterns, denoted "DeBruijn(Patterns)",
by gluing identically labeled nodes in CompositionGraph(Patterns),
which yields the following algorithm.

    DEBRUIJN(Patterns)
        represent every k-mer in Patterns as an isolated edge between its prefix and suffix
        glue all nodes with identical labels, yielding the graph DeBruijn(Patterns)
        return DeBruijn(Patterns)

"De Bruijn Graph" from k-mers Problem
Construct the "de Bruijn graph" from a collection of k-mers.

Given: A collection of k-mers Patterns.

Return: The "de Bruijn graph" "DeBruijn(Patterns)", in the form of an adjacency list.

Sample Dataset
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG

Sample Output
AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> PREV: Construct De Bruijn Graph (BA3D)
      -> HERE: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> NEXT: Eulerian Cycle (BA3F)

═════════════════════════════════════════════════

References:
- How to find if a given key exists in a C++ std::map
  https://stackoverflow.com/questions/1939953/how-to-find-if-a-given-key-exists-in-a-c-stdmap
*/

#include <unordered_set>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

// BA03E: De Bruijn graph from k-mer patterns
unordered_map<string, vector<string>> debruijn_from_kmers(vector<string> const& patterns) {
    unordered_map<string, vector<string>> map;
    for (auto pattern : patterns) {
        string prefix = pattern.substr(0, pattern.length() - 1);
        string suffix = pattern.substr(1);
        if (map.find(prefix) == map.end()) {  // if key is not in map
            map.insert({ prefix, {suffix} });
        }
        else { // if key is in map
            map.at(prefix).push_back(suffix);
        }
    }
    return map;
}

// BA03E: De Bruijn graph from k-mer patterns, using set  <--- Don't use this fx. This is wrong. TODO: Why?
unordered_map<string, unordered_set<string>> debruijn_from_kmers_set(vector<string> const& patterns) {
    unordered_map<string, unordered_set<string>> map;
    for (auto pattern : patterns) {
        string prefix = pattern.substr(0, pattern.length() - 1);
        string suffix = pattern.substr(1);
        if (map.find(prefix) == map.end()) {  // if key is not in map
            map.insert({ prefix, {suffix} });
        }
        else { // if key is in map
            map.at(prefix).insert(suffix);
        }
    }
    return map;
}

// main
int main() {
    vector<string> const lines = read_lines("/home/wsl/rosalind/data/ba03e.txt");

    // print to screen
    clock_t tStart = clock();
    print_unordered_map(debruijn_from_kmers(lines));
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    /*
    // print result to text file
    ofstream ofs{ "/home/wsl/rosalind/data/ba03e_output.txt" };
    auto cout_buff = cout.rdbuf();
    cout.rdbuf(ofs.rdbuf());
    print_unordered_map(debruijn_from_kmers(lines));
    cout.rdbuf(cout_buff);
    */
}