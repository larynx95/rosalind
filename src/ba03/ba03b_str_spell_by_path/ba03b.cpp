/*
Rosalind: BA3B
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern[1], ... , Pattern[n]
such that the last k - 1 symbols of Pattern[i] are equal
to the first k - 1 symbols of Pattern[i+1] for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> PREV: Generate the k-mer Composition of a String (BA3A)
      -> HERE: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * NEXT: Construct OverapGraph (BA3C)

References:
- Best way to extract a subvector from a vector?
  https://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

// BA03B: string spell by ordered genome path
string str_spell(vector<string> patterns) {
    string result = patterns[0];
    vector<string> subvec{ patterns.begin() + 1, patterns.end() };
    for (auto& pattern : subvec) {
        result += pattern.substr(pattern.length() - 1, 1);
    }
    return result;
}

// BA03B: string spell by ordered genome path
string str_spell2(vector<string> patterns) {
    string result = patterns[0];
    vector<string> subvec{ patterns.begin() + 1, patterns.end() };
    for (auto it = patterns.begin() + 1; it != patterns.end(); it++) {
        result += (*it).substr(it->length() - 1, 1);  // two different ways of dereferencing
    }
    return result;
}

// main
int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03b.txt");

    clock_t tStart = clock();
    cout << str_spell2(lines) << endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}
