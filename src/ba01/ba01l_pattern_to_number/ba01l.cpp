/*
Rosalind: BA1L
Implement PatternToNumber

Implement PatternToNumber
Convert a DNA string to a number.

Given: A DNA string Pattern.
Return: PatternToNumber(Pattern).

Sample Dataset
AGT

Sample Output
11

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
- modulo operator
- number types: long long int, long double
- std::pow(), return type (Be careful)
*/


#include <cmath>
#include "/home/wsl/rosalind/include/utils.hpp"

// BA01L: pattern to number
long long int pattern_to_number(std::string pattern) {
    long long int number = 0;
    int len = pattern.size();
    for (int i = 0; i < len; i++) {
        int pw = len - i - 1;
        int nuc;
        switch (pattern[i]) {
        case 'A':
            nuc = 0;
            break;
        case 'C':
            nuc = 1;
            break;
        case 'G':
            nuc = 2;
            break;
        case 'T':
            nuc = 3;
            break;
        }
        long double pw_val = std::pow(4.0, pw * 1.0);
        number += (pw_val * nuc);
    }
    return number;
}

// main
int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01l.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];

    std::cout << pattern_to_number(pattern) << std::endl;

    return 0;
}
