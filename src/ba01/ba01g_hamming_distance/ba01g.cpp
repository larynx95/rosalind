/*
Rosalind: BA1G
Compute the Hamming Distance Between Two Strings

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
For example, CGAAT and CGGAC have two mismatches. The number of mismatches
between strings p and q is called the Hamming distance between these strings and
is denoted HammingDistance(p, q).

Hamming Distance Problem
Compute the Hamming distance between two DNA strings.

Given: Two DNA strings.

Return: An integer value representing the Hamming distance.

Sample Dataset
GGGCCGTTGGT
GGACCGTTGAC

Sample Output
3

═════════════════════════════════════════════════

Plan 1.
- simple loop

Topic:
- Hamming distance
- TODO: Better alogithms?
*/

#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

template <typename T>
bool is_vec_contains(std::vector<T> vec, T elem);
int get_hamming_distance(std::string const& str1, std::string const& str2);  // not change, so use 'const&'

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01g.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string str1 = lines[0];
    std::string str2 = lines[1];

    clock_t tStart3 = clock();
    int hamming_distance = get_hamming_distance(str1, str2);
    std::cout << "Hamming distance: " << hamming_distance << std::endl;
    printf("Time taken: %.6fs\n", (double)(clock() - tStart3) / CLOCKS_PER_SEC);

    return 0;
}

/************************************************
definition
************************************************/

// BA01G: hamming distance
int get_hamming_distance(std::string const& str1, std::string const& str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}
