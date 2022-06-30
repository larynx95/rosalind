/*
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif
if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

╔════════════════════════════════════════════════════════════════════════════════════╗
║ MOTIFENUMERATION(Dna, k, d)                                                        ║
║     Patterns <- an empty set                                                       ║
║     for each k-mer Pattern in Dna                                                  ║
║         for each k-mer Pattern' differing from Pattern by at most d mismatches     ║
║             if Pattern' appears in each string from Dna with at most d mismatches  ║
║                 add Pattern' to Patterns                                           ║
║     remove duplicates from Patterns                                                ║
║     return Patterns                                                                ║
╚════════════════════════════════════════════════════════════════════════════════════╝

Implanted Motif Problem
Implement MotifEnumeration (shown above) to find all (k, d)-motifs
in a collection of strings.

Given: Integers k and d, followed by a collection of strings Dna.
Return: All (k, d)-motifs in Dna.

Sample Dataset
3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output
ATA ATT GTT TTT

═════════════════════════════════════════════════

Info.
  * using BA1N neighbors function, set intersection
    a. get all d-neighbors from all k-mers of the first DNA string
    b. repeat set intersection between the first d-neighbors and each of d-neighbor from each DNA string
  * sample dataset
    dnas     k=3    d=1                                      motifs
    -----------------------------------------------------------------
    ATTTGGC  ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT  ATA ATT GTT TTT
             TTT    TTA TTC TTG TAT TCT TGT ATT CTT GTT TTT
             TTG    TTA TTC TAG TCG TGG ATG CTG GTG TTG TTT
             TGG    TGA TGC TAG TCG AGG CGG GGG TGG TTG TGT
             GGC    GGA GAC GCC AGC CGC GGC TGC GTC GGG GGT
    TGCCTTA  TGC    TGA TAC TCC AGC CGC GGC TGC TTC TGG TGT
             GCC    GCA GAC ACC CCC GCC TCC GGC GTC GCG GCT
             CCT    CCA CCC CCG CAT ACT CCT GCT TCT CGT CTT
             CTT    CTA CTC CTG CAT CCT CGT ATT CTT GTT TTT
             TTA    TAA TCA TGA ATA CTA GTA TTA TTC TTG TTT
    CGGTATC  CGG    CGA CGC CAG CCG AGG CGG GGG TGG CTG CGT
             GGT    GGA GGC GGG GAT GCT AGT CGT GGT TGT GTT
             GTA    GAA GCA GGA ATA CTA GTA TTA GTC GTG GTT
             TAT    TAA TAC TAG AAT CAT GAT TAT TCT TGT TTT
             ATC    ATA AAC ACC AGC ATC CTC GTC TTC ATG ATT
    GAAAATT  GAA    AAA CAA GAA TAA GCA GGA GTA GAC GAG GAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAT    AAA AAC AAG AAT CAT GAT TAT ACT AGT ATT
             ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT

═════════════════════════════════════════════════

References:
- Slicing a vector in C++
  https://stackoverflow.com/questions/50549611/slicing-a-vector-in-c
- How to find the intersection of two STL sets?
  https://stackoverflow.com/questions/13448064/how-to-find-the-intersection-of-two-stl-sets
- Best way to extract a subvector from a vector?
  https://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
- TODO: set intersection, modifying iterator during iteration
  - How to check that an element is in a set?
    https://stackoverflow.com/questions/1701067/how-to-check-that-an-element-is-in-a-stdset
  - delete elem in a set: 'set.erase'
  - ☆☆☆☆☆ How to remove elements from an std::set while iterating over it
    https://stackoverflow.com/questions/20627458/how-to-remove-elements-from-an-stdset-while-iterating-over-it
*/

#include <unordered_set>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

int hamming_distance(string const& str1, string const& str2);
vector<string> neighbors_rec(string, int);
unordered_set<string> motif_enumeration(vector<string> dnas, int k, int d);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const string filename = argc < 2 ? "/home/wsl/rosalind/data/ba02a.txt" : argv[1];
    vector<string> lines = read_lines(filename);
    vector<string> splitted = split_str(lines[0], " ");
    int k = stoi(splitted[0]);
    int d = stoi(splitted[1]);
    vector<string> dnas;
    dnas.assign(lines.begin() + 1, lines.end());

    clock_t tStart = clock();
    for (auto motif : motif_enumeration(dnas, k, d)) {
        cout << motif << " ";
    }
    cout << endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // practice
    // modifying original iterator during iteration
    unordered_set<string> a = { "AAT", "ACT", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CGC", "CGG", "CTG", "CTT", "GAC", "GCC", "GGA", "GGC", "GGG", "GGT", "GTC", "GTG", "GTT", "TAG", "TAT", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT" };
    unordered_set<string> b = { "ACC","ACT","AGC","ATA","ATT","CAT","CCA","CCC","CCG","CCT","CGC","CGT","CTA","CTC","CTG","CTT","GAC","GCA","GCC","GCG","GCT","GGC","GTA","GTC","GTT","TAA","TAC","TCA","TCC","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT" };
    unordered_set<string> c = { "AAC","AAT","ACC","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAG","CAT","CCG","CGA","CGC","CGG","CGT","CTA","CTC","CTG","GAA","GAT","GCA","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCT","TGG","TGT","TTA","TTC","TTT" };
    unordered_set<string> e = { "AAA","AAC","AAG","AAT","ACA","ACT","AGA","AGT","ATA","ATC","ATG","ATT","CAA","CAT","CTT","GAA","GAC","GAG","GAT","GCA","GGA","GTA","GTT","TAA","TAT","TTT" };

    for (auto it = a.begin(); it != a.end(); ) {
        if (b.find(*it) == b.end()) {
            it = a.erase(it);
        }
        else {
            it++;
        }
    }

    for (auto it = a.begin(); it != a.end(); ) {
        if (c.find(*it) == c.end()) {
            it = a.erase(it);
        }
        else {
            it++;
        }
    }

    for (auto it = a.begin(); it != a.end(); ) {
        if (e.find(*it) == e.end()) {
            it = a.erase(it);
        }
        else {
            it++;
        }
    }

    for (auto e : a) {
        cout << e << " ";
    }
    cout << endl;



    return 0;
}

/************************************************
definitions
************************************************/

// BA01G: hamming distance
int hamming_distance(string const& str1, string const& str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

// BA01N: recursive neighbors, return unordered_set
unordered_set<string> neighbors(string pattern, int d) {
    if (d == 0) {
        return { pattern };
    }
    if (pattern.length() == 1) {
        return { "A", "C", "G", "T" };
    }
    unordered_set<string> neighborhood;
    unordered_set<string> suffix_neighbors = neighbors(pattern.substr(1), d);
    for (auto text : suffix_neighbors) {
        if (hamming_distance(pattern.substr(1), text) < d) {
            for (auto nuc : { "A", "C", "G", "T" }) {  // "ACGT" does not work, TODO: Why?
                neighborhood.insert(nuc + text);
            }
        }
        else {
            neighborhood.insert(pattern[0] + text);
        }
    }
    return neighborhood;
}

// BA02A: motif enumeration
unordered_set<string> motif_enumeration(vector<string> dnas, int k, int d) {
    unordered_set<string> candidates;
    // get all d-neighbors from the first dna string: candidates
    for (int i = 0; i < dnas[0].length() - k + 1; i++) {
        unordered_set<string> nbs = neighbors(dnas[0].substr(i, k), d);
        candidates.insert(nbs.begin(), nbs.end());
    }
    // get all d-neighbors from each dna string, compare, do intersection
    for (int i = 1; i < dnas.size(); i++) {
        // get all d-neighbors from i-th dna string
        string dna = dnas[i];
        unordered_set<string> ith_nbs;
        for (int j = 0; j < dna.length() - k + 1;j++) {
            unordered_set<string> nbs = neighbors(dna.substr(j, k), d);
            ith_nbs.insert(nbs.begin(), nbs.end());
        }
        // unordered_set intersection, TODO: I wasted hours fixing this part.
        // ★★★★★ iterator changes during iteration problem!
        for (auto it = candidates.cbegin(); it != candidates.cend();) {
            if (ith_nbs.find(*it) == ith_nbs.end()) {
                it = candidates.erase(it);
            }
            else {
                it++;
            }
        }
    }
    return candidates;
}

