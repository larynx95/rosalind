/*
Rosalind: BA1D
Find All Occurrences of a Pattern in a String

In this problem, we ask a simple question: how many times can one string occur
as a substring of another? Recall from “Find the Most Frequent Words in a
String” that different occurrences of a substring can overlap with each other.
For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
Find all occurrences of a pattern in a string.

Given: Strings Pattern and Genome.

Return: All starting positions in Genome where Pattern appears as a substring.
Use 0-based indexing.

Sample Dataset
ATAT
GATATATGCATATACTT

Sample Output
1 3 9

═════════════════════════════════════════════════

*/

#include "/home/wsl/rosalind/include/utils.hpp"

// BA01D: find all occurrences
std::vector<int> find_all_ocurrences(std::string pattern, std::string genome) {
    std::vector<int> indices;
    int len_genome = genome.length();
    int len_pattern = pattern.length();

    for (int i = 0; i < len_genome - len_pattern + 1; i++) {
        if (genome.substr(i, len_pattern) == pattern) {
            indices.push_back(i);
        }
    }
    return indices;
}

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01d.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    std::string genome = lines[1];

    std::vector<int> indices = find_all_ocurrences(pattern, genome);
    print_vec(indices);

    return 0;
}
