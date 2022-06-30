/*
Rosalind: BA1M
Implement NumberToPattern

Convert an integer to its corresponding DNA string.

Given: Integers index and k.
Return: NumberToPattern(index, k).

Sample Dataset
45
4

Sample Output
AGTC

═════════════════════════════════════════════════
Pseudocode:

  ╔════════════════════════════════════════════════════════════════════╗
  ║  PATTERNTONUMBER(Pattern)                                          ║
  ║      if Pattern contains no symbols                                ║
  ║          return 0                                                  ║
  ║      symbol LASTSYMBOL(Pattern)                                    ║
  ║      Prefix PREFIX(Pattern)                                        ║
  ║      return 4 · PATTERNTONUMBER(Prefix) + SYMBOLTONUMBER(symbol)   ║
  ╚════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NUMBERTOPATTERN(index , k)                              ║
  ║      if k = 1                                            ║
  ║          return NUMBERTOSYMBOL(index)                    ║
  ║      prefixIndex QUOTIENT(index, 4)                      ║
  ║      r REMAINDER(index, 4)                               ║
  ║      symbol NUMBERTOSYMBOL(r)                            ║
  ║      PrefixPattern NUMBERTOPATTERN(prefixIndex, k  1)    ║
  ║      return concatenation of PrefixPattern with symbol   ║
  ╚══════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

Plan 1.
- find a rule between given string and index
kmer    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
index    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
x % 4    0  1  2  3  0  1  2  3  0  1  2  3  0  1  2  3
x / 4    0  0  0  0  1  1  1  1  2  2  2  2  3  3  3  3

              A:0, C:1, G:2, T:3

  exp   1 0              2 1 0
  Nuc   T T : 3 3        G A T : 2 0 3
        │ └── 4^0*3      │ │ └── 4^0*2
        └──── 4^1*3      │ └──── 4^1*0
               = 15      └────── 4^2*3
                                  = 50

Topic:
- '/', '%' operators
- TODO: Implement recursive number_to_pattern fx.
*/

#include <cmath>
#include "/home/wsl/rosalind/include/utils.hpp"

// BA01M: number to pattern
std::string number_to_pattern(int number, int num_k) {
    std::string pattern;
    for (int i = 0; i < num_k; i++) {
        int temp = number / std::pow(4, num_k - i - 1);
        switch (temp) {
        case 0:
            pattern += 'A';
            break;
        case 1:
            pattern += 'C';
            break;
        case 2:
            pattern += 'G';
            break;
        case 3:
            pattern += 'T';
            break;
        }
        number -= temp * std::pow(4, num_k - i - 1);
    }
    return pattern;
}

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01m.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    int number = std::stoi(lines[0]);
    int num_k = std::stoi(lines[1]);

    std::cout << number_to_pattern(number, num_k) << std::endl;

    return 0;
}