/*
Rosalind: BA3A (difficulty: 1/5)
Generate the k-mer Composition of a String

Given a string Text, its k-mer composition Compositionk(Text)
is the collection of all k-mer substrings of Text
(including repeated k-mers).
For example,

Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}

Note that we have listed k-mers in lexicographic order
(i.e., how they would appear in a dictionary)
rather than in the order of their appearance in TATGGGGTGC.
We have done this because the correct ordering of the reads
is unknown when they are generated.

String Composition Problem
Generate the k-mer composition of a string.

Given: An integer k and a string Text.

Return: Compositionk(Text)
        (the k-mers can be provided in any order).

Sample Dataset
5
CAATCCAAC

Sample Output
AATCC
ATCCA
CAATC
CCAAC
TCCAA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> HERE: Generate the k-mer Composition of a String (BA3A)
      -> NEXT: String Reconstruction Problem with k-mer composition(BA3B)
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

vector<string> all_kmers(string const& text, int const k);

/************************************************
main
************************************************/

int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03a.txt");
    int k = stoi(lines[0]);
    string text = lines[1];

    clock_t tStart = clock();
    for (auto kmer : all_kmers(text, k)) {
        cout << kmer << endl;
    }
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

/************************************************
definitions
************************************************/

// BA03A: get all k-mers from a text
vector<string> all_kmers(string const& text, int const k) {
    vector<string> kmers;
    for (size_t i = 0; i < text.length() - k + 1; i++) {
        kmers.push_back(text.substr(i, k));
    }
    return kmers;
}
