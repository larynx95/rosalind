/*
Rosalind: BA1N
Generate the d-Neighborhood of a String

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming
distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String
Find all the neighbors of a pattern.

Given: A DNA string Pattern and an integer d.
Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
ACG
1

Sample Output
CCG
TCG
GCG
AAG
ATG
AGG
ACA
ACC
ACT
ACG

═════════════════════════════════════════════════
Pseudocode:
  * immediate neighbors
    ACG, d=1
    A --> CCG, GCG, TCG
    C --> AAG, AGG, ATG
    G --> ACA, ACC, ACT
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
  A C G         original pattern: ACG
  │ │ └─ substitute with A, C, T: ACA, ACC, ACT
  │ └─── substitute with A, G, T: AAG, AGG, ATG
  └───── substitute with C, G, T: CCG, GCG, TCG
  result: ACG, ACA, ACC, ACT, AAG, AGG, ATG, CCG, GCG, TCG

  * iterative neighbors
    ╔══════════════════════════════════════════════════════════════╗
    ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
    ║    Neighborhood <- set consisting of single string Pattern   ║
    ║    for j = 1 to d                                            ║
    ║      for each string Pattern' in Neighborhood                ║
    ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
    ║        remove duplicates from Neighborhood                   ║
    ║    return Neighborhood                                       ║
    ╚══════════════════════════════════════════════════════════════╝

  * neighbors - recursive function
    (1) base case
      - if d == 0          => {pattern}
      - if |pattern| == 0  => {}
      - if |pattern| == 1  => {'A', 'C', 'G', 'T'}
    (2) recursive case
      - if neighbors(pattern[1:]) < d   => can prepend any nucleotide
      - if neighbors(pattern[1:]) >= d  => can change the first nucleotide
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

Info.
  [d-1] Neighbors: get d-1 neighbers from original pattern
    ACG ┌ [A]CG : [C]CG, [G]CG, [T]CG    --- (1)
        │
        │
        │ A[C]G : A[A]G, A[G]G, A[T]G    --- (2)
        │
        │
        └ AC[G] : AC[A], AC[C], AC[T]    --- (3)

  [d-2] Neighbors: get d-1 neighbors from previous step, then concetanate
    In case (1):
      - substring is "CG" (pattern').
      - get (d-1 neighbors)' from "CG": AG, GG, TG, CA, CC, CT
      - concetenate char + each of (d-1 neighbors)'
        ex) [C]CG:  CAG, CGG, CTG, CCA, CCC, CCT
    In case (2):
      - substring is splitted by the character in previous d-1 step
      ┌ Additional step:
      │ - First, we should combine prefix and postfix togerther.
      │   ex) A[A]G: "A" + "G" = "AG"
      │ - get (d-1 neighbors)' from "AG": CG, GG, TG, AA, AC, AT
      │ Option A: concat -> split -> append all):
      │   - remember: index of middle char, the sizes of pre/postfix
      │   - split each (d-1 neighbor)' into prefix and postfix parts
      │   - concatenate prefix + middle char + postfix for each (d-1 neighbor)
      │     ex) A[A]G: C+A+G, G+A+G, T+A+G, A+A+A, A+A+C, A+A+T
      │                CAG    GAG    TAG    AAA    AAC    AAT
      │ Option B. concat -> insert
      │   - remember: index of middle char
      │   - insert middle char into each of (d-1 neighbor)' at given index
      │     ex) A[A]G: C^G, G^G, T^G, A^A, A^C, A^T
      └                CAG, GAG, TAG, AAA, AAC, AAT
    In case (3):
      - Actually, (concat -> split -> append) maybe equal to (concat -> insert)
      - Actually, in all three cases, you may have to go through the same process.

═════════════════════════════════════════════════
References:
- Parse (split) a string in C++ using string delimiter (standard C++)
  https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
- Concatenating two std::vectors
  https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
  vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
- How to convert vector to set?
  https://stackoverflow.com/questions/20052674/how-to-convert-vector-to-set
  std::vector<std::string> v;
  std::set<std::string> s(v.begin(), v.end());
- C++, copy set to vector
  https://stackoverflow.com/questions/5034211/c-copy-set-to-vector
- What's the most efficient way to erase duplicates and sort a vector?
  https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
- C++ equivalent of Python String Slice?
  https://stackoverflow.com/questions/27992264/c-equivalent-of-python-string-slice
- "ACGT" vs {"A", "C", "G", "T"}
- vector literal; {"A", "C", "G", "T"}
- Benefits of vector<char> over string?
  https://stackoverflow.com/questions/11358879/benefits-of-vectorchar-over-string
- append set to another set
  https://stackoverflow.com/questions/2607235/append-set-to-another-set
- ♨ modifying iterator during iteration
  - How to properly iterate over container while modifying it? [duplicate]
    https://stackoverflow.com/questions/59667104/how-to-properly-iterate-over-container-while-modifying-it
  - erase set element while iterating///
    https://stackoverflow.com/questions/51208537/erase-set-element-while-iterating
  - How to remove from a map while iterating it?
    https://stackoverflow.com/questions/8234779/how-to-remove-from-a-map-while-iterating-it/8234813#8234813
*/

#include <cmath>
#include <unordered_set>  // Don't use <set> here! unexpected result!
#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
function prototypes
************************************************/

template <typename T>
void print_vec_one_line(const T& vec);
int hamming_distance(std::string const& str1, std::string const& str2);

std::vector<std::string> neighbors_immediate(std::string const&);
std::vector<std::string> neighbors_iter(std::string const&, int const);
std::vector<std::string> neighbors_rec(std::string, int);

std::unordered_set<std::string> neighbors_immediate_set(std::string);
std::unordered_set<std::string> neighbors_iter_set(std::string, int);
std::unordered_set<std::string> neighbors_rec_set(std::string, int);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01n.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    int d = std::stoi(lines[1]);
    d = 2;

    clock_t tStart = clock();
    std::cout << "--- neighbors, immediate --- " << std::endl;
    std::vector<std::string> imm = neighbors_immediate(pattern);
    print_vec_one_line(imm);
    std::cout << "size: " << imm.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    tStart = clock();
    std::cout << "\n------- neighbors, immediate, set ------- " << std::endl;
    std::unordered_set<std::string> nb_imm_set = neighbors_immediate_set(pattern);
    print_vec_one_line(nb_imm_set);
    std::cout << "size: " << nb_imm_set.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    /*
    tStart = clock();
    std::cout << "\n------- neighbors, iterative ------- " << std::endl;
    std::vector<std::string> nb_iter = neighbors_iter(pattern, d);
    print_vec_one_line(nb_iter);
    std::cout << "size: " << nb_iter.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    */

    /*
    tStart = clock();
    std::cout << "\n------- neighbors, iterative, set ------- " << std::endl;
    std::unordered_set<std::string> nb_iter_set = neighbors_iter_set(pattern, d);
    print_vec_one_line(nb_iter_set);
    std::cout << "size: " << nb_iter_set.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    */

    tStart = clock();
    std::cout << "\n------- neighbors, recursive ------- " << std::endl;
    std::vector<std::string> nb_rec = neighbors_rec(pattern, d);
    print_vec_one_line(nb_rec);
    std::cout << "size: " << nb_rec.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    tStart = clock();
    std::cout << "\n------- neighbors, recursive, set ------- " << std::endl;
    std::unordered_set<std::string> nb_rec_set = neighbors_rec_set(pattern, d);
    print_vec_one_line(nb_rec_set);
    std::cout << "size: " << nb_rec_set.size() << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}

/************************************************
definitions
************************************************/

// helper fx.: print vector
template <typename T>
void print_vec_one_line(const T& vec) {
    for (auto iter = vec.begin(); iter != vec.end(); iter++) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}

// BA01G: hamming distance
int hamming_distance(std::string const& str1, std::string const& str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

// BA01N: immediate neighbors
std::vector<std::string> neighbors_immediate(std::string const& pattern) {
    std::vector<std::string> vec;
    for (int i = 0; i < pattern.size(); i++) {
        char symbol = pattern[i];
        std::string nucleotides = "ACGT";
        for (auto nuc : nucleotides) {
            if (symbol != nuc) {
                std::string pre = pattern.substr(0, i);
                std::string post = pattern.substr(i + 1, pattern.size() - i - 1);
                std::string temp = pre + nuc + post;
                vec.push_back(temp);
            }
        }
    }
    vec.push_back(pattern);
    return vec;
}

// BA01N: neighbors, immediate, returning Set
std::unordered_set<std::string> neighbors_immediate_set(std::string pattern) {
    std::unordered_set<std::string> neighbors{ pattern };
    for (int i = 0; i < pattern.size(); i++) {
        char symbol = pattern[i];
        std::string nucleotides = "ACGT";
        for (auto chr : nucleotides) {
            if (symbol != chr) {
                std::string pre = pattern.substr(0, i);
                std::string post = pattern.substr(i + 1, pattern.size() - i - 1);
                std::string temp = pre + chr + post;
                neighbors.insert(temp);
            }
        }
    }
    return neighbors;
}

/*
// BA01N: neighbors, iterative TODO: Fix this!
std::vector<std::string> neighbors_iter(std::string const& pattern, int const d) {
    std::vector<std::string> neighbors{ pattern };
    for (int i = 0; i < d; i++) {
        for (auto pat : neighbors) {  // modifying iterator during iteration, TODO: Is this right?
            // for (auto it = neighbors.cbegin(); it != neighbors.cend();) {
                // get immediate neighbors
            std::vector<std::string> imm = neighbors_immediate(pat);
            // add immediate neighbors to neigbors
            neighbors.insert(neighbors.end(), imm.begin(), imm.end());
            // remove duplicates using set
            std::unordered_set<std::string> set_temp;
            for (unsigned i = 0; i < neighbors.size(); i++) set_temp.insert(neighbors[i]);
            neighbors.assign(set_temp.begin(), set_temp.end());
        }
    }
    return neighbors;
}
*/

/*
// BA01N: neighbors, iterative, returning Set, TODO: fix this
std::unordered_set<std::string> neighbors_iter_set(std::string pattern, int d) {
    std::unordered_set<std::string> neighbors{ pattern };
    for (int i = 0; i < d; i++) {
        // for (auto pat : neighbors) {  // unexpected result
        for (auto it = neighbors.cbegin(); it != neighbors.cend();) {
            std::unordered_set<std::string> imm = neighbors_immediate_set(*it);// get immediate neighbors
            neighbors.insert(imm.begin(), imm.end());// add immediate neighbors to neigbors
        }
    }
    return neighbors;
}
*/

// BA01N: neighbors, recursive
std::vector<std::string> neighbors_rec(std::string pattern, int d) {
    if (d == 0) {
        return { pattern };
    }
    if (pattern.length() == 1) {
        return { "A", "C", "G", "T" };
    }
    std::unordered_set<std::string> neighborhood;
    std::vector<std::string> suffix_neighbors = neighbors_rec(pattern.substr(1), d);
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
    std::vector<std::string> vec(neighborhood.begin(), neighborhood.end());
    return vec;
}

// BA01N: neighbors, recursive, returning Set
std::unordered_set<std::string> neighbors_rec_set(std::string pattern, int d) {
    if (d == 0) {
        return { pattern };
    }
    if (pattern.length() == 1) {
        return { "A", "C", "G", "T" };
    }
    std::unordered_set<std::string> neighborhood;
    std::unordered_set<std::string> suffix_neighbors = neighbors_rec_set(pattern.substr(1), d);
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