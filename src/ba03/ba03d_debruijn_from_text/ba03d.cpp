/*
Rosalind: BA3D
Construct the "De Bruijn[broin]" Graph of a String

Given a genome Text, PathGraph_k(Text) is the path consisting of |Text| - k + 1 edges,
where the i-th edge of this path is labeled by the i-th k-mer in Text
and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text.
The de Bruijn graph DeBruijn_k(Text) is formed
by gluing identically labeled nodes in PathGraph_k(Text).

De Bruijn Graph from a String Problem
Construct the de Bruijn graph of a string.

Given: An integer k and a string Text.

Return:DeBruijn_k(Text), in the form of an adjacency list.

Sample Dataset
4
AAGATTCTCTAC

Sample Output (1:many)
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC
TTC -> TCT

═════════════════════════════════════════════════

References:
- Check if a given key exists in a map or not in C++
  https://www.techiedelight.com/check-given-key-exists-map-cpp/
- What is the quickest way of inserting/updating std::unordered_map elements without using an if?
  https://stackoverflow.com/questions/19197799/what-is-the-quickest-way-of-inserting-updating-stdunordered-map-elements-witho
- std::unordered_map<Key,T,Hash,KeyEqual,Allocator>::insert
  https://en.cppreference.com/w/cpp/container/unordered_map/insert
- Traversing a map (or unordered_map) in C++ STL
  https://www.geeksforgeeks.org/traversing-a-map-or-unordered_map-in-cpp-stl/
- Iterating over unordered_map of vectors
  https://stackoverflow.com/questions/62614386/iterating-over-unordered-map-of-vectors
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

// BA03D: De Bruijn graph from a text
unordered_map<string, vector<string>> debruijn_from_text(string text, int k) {
    unordered_map<string, vector<string>> graph;
    for (size_t i = 0; i < text.length() - k + 1; i++) {
        string prefix = text.substr(i, k - 1);
        string suffix = text.substr(i + 1, k - 1);
        if (graph.find(prefix) == graph.end()) { // if key is not in unordered_map
            graph.insert({ prefix, {suffix} });
        }
        else { // if key is in unordered_map
            graph.at(prefix).push_back(suffix);
        }
    }
    return graph;
}

// main
int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03d.txt");
    int k = stoi(lines[0]);
    string text = lines[1];

    // print result to screen
    clock_t tStart = clock();
    for (const auto& [key, values] : debruijn_from_text(text, k)) {
        cout << key << " -> ";
        string val;
        for (auto it = values.begin(); it != values.end(); it++) {
            if (it != values.begin()) {
                val += ",";
            }
            val += *it;
        }
        cout << val << endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    /*
    // print result to text file
    ofstream ofs{ "/home/wsl/rosalind/data/ba03d_output.txt" };
    auto cout_buff = cout.rdbuf();
    cout.rdbuf(ofs.rdbuf());
    for (const auto& [key, values] : debruijn_from_text(text, k)) {
        cout << key << " -> ";
        string val;
        for (auto it = values.begin(); it != values.end(); it++) {
            if (it != values.begin()) {
                val += ",";
            }
            val += *it;
        }
        cout << val << endl;
    }
    cout.rdbuf(cout_buff);
    */
}