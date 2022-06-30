/*
Rosalind: BA1F
Find a Position in a Genome Minimizing the Skew

Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference
between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi
(Genome) denote the prefix (i.e., initial substring) of Genome of length i. For
example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i
(from 0 to |Genome|).

Sample Dataset
CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

Sample Output
53 97

═════════════════════════════════════════════════

Plan 1.
- get GC skew table, find min value
- Find the indices of the values ​​that match the min value

Topics:
- GC skewness
*/

#include "/home/wsl/rosalind/include/utils.hpp"

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01f.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string genome = lines[0];
    std::vector<int> skew_table = { 0 };
    int skew = 0;
    int min_skew = 0;
    for (auto nucleotide : genome) {
        if (nucleotide == 'C') {
            skew -= 1;
        }
        if (nucleotide == 'G') {
            skew += 1;
        }
        skew_table.push_back(skew);
        if (skew < min_skew) min_skew = skew;
    }
    std::vector<int> indices;
    for (int i = 0; i < skew_table.size(); i++) {
        if (skew_table[i] == min_skew) indices.push_back(i);
    }
    print_vec(indices);
    return 0;
}