/*
Rosalind: BA1A
Compute the Number of Times a Pattern Appears in a Text

This is the first problem in a collection of "code challenges"
to accompany Bioinformatics Algorithms:
An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

A k-mer is a string of length k.
We define Count(Text, Pattern) as the number of times
that a k-mer Pattern appears as a substring of Text.
For example,

Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.

We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2)
since we should account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern),
our plan is to “slide a window” down Text,
checking whether each k-mer substring of Text matches Pattern.
We will therefore refer to the k-mer starting at position i of Text as Text(i,
k). Throughout this book, we will often use 0-based indexing, meaning that we
count starting at 0 instead of 1. In this case, Text begins at position 0 and
ends at position |Text| − 1 (|Text| denotes the number of symbols in Text). For
example, if Text = GACCATACTG, then Text(4, 3) = ATA. Note that the last k-mer
of Text begins at position |Text| − k, e.g., the last 3-mer of GACCATACTG starts
at position 10 −3 = 7. This discussion results in the following pseudocode for
computing Count(Text, Pattern).

╔═════════════════════════════════════════╗
║ PatternCount(Text, Pattern)             ║
║     count <- 0                          ║
║     for i <- 0 to |Text| - |Pattern|    ║
║         if Text(i, |Pattern|) = Pattern ║
║             count <- count + 1          ║
║     return count                        ║
╚═════════════════════════════════════════╝

Implement PatternCount
Given: {DNA strings}} Text and Pattern.
Return: Count(Text, Pattern).

Sample Dataset
GCGCG
GCG

Sample Output
2

═════════════════════════════════════════════════
References:
- Read file line by line using ifstream in C++
  https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c
- counting the number of lines in a text file
  https://stackoverflow.com/questions/3482064/counting-the-number-of-lines-in-a-text-file
- C++ String Array, Loading lines of text from file
  https://stackoverflow.com/questions/3901977/c-string-array-loading-lines-of-text-from-file
*/

#include <chrono>
#include "/home/wsl/rosalind/include/utils.hpp"

// BA01A: pattern count
int pattern_count(std::string text, std::string pattern) {
    int count = 0;
    int k = pattern.length();
    for (int i = 0; i < text.length() - k + 1; i++) {
        std::string kmer = text.substr(i, k);
        if (kmer == pattern) count++;
    }
    return count;
}

// main
int main(int argc, const char* argv[]) {
    // if argc == 2, filename is argv[1]
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01a.txt" : argv[1];
    // error handling, if no such file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cout << "no file" << std::endl;
    }
    // read text file, store data into array of strings
    std::vector<std::string> lines = read_lines(filename);
    // given data: 'pattern' and 'text'
    std::string text{ lines[0] };
    std::string pattern{ lines[1] };

    // check time
    clock_t tStart = clock();
    printf("%d\n", pattern_count(text, pattern));
    printf("Time taken: %.5fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}
