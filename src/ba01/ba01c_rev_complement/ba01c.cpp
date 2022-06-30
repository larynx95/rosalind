/*
Rosalind: BA1C
Find the Reverse Complement of a String

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C'
and 'G'. Given a nucleotide p, we denote its complementary nucleotide as p. The
reverse complement of a DNA string Pattern = p1…pn is the string rPattern = pn …
p1 formed by taking the complement of each nucleotide in Pattern, then reversing
the resulting string.

For example, the reverse complement of Pattern = "GTCA" is rPattern = "TGAC".

Reverse Complement Problem
Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: rPattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT

═════════════════════════════════════════════════

Plan 1.
  - reverse string, then replace each nucleotide with its couterpart

Topics:
  - change character in string
  - call by value, call by refererence
  - swap
  - reverse array, vector

═════════════════════════════════════════════════

*/

#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"

// BA01C: reverse complement
std::string rev_complement(std::string& pattern) {
    std::string complement;
    for (auto rev_iter = pattern.rbegin(); rev_iter != pattern.rend();  // reverse iterator
        rev_iter++) {
        char temp = *rev_iter;
        switch (temp) {
        case 'A':
            complement += 'T';
            break;
        case 'C':
            complement += 'G';
            break;
        case 'G':
            complement += 'C';
            break;
        case 'T':
            complement += 'A';
            break;
        }
    }
    return complement;
}

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01c.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string text = lines[0];

    clock_t tStart = clock();
    std::cout << rev_complement(text) << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}


