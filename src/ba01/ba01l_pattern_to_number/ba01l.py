"""
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

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> find frequent words in a string (BA1B)
      -> find all occurrence of a pattern in a string (BA1D)
      -> Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> hamming distance (BA1G)
      -> PREV: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> NEXT: frequency array
          -> HERE: number to pattern (BA1M), pattern to number (BA1L)

Plan 1.
- pseudocode:
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

"""

import time


# BA01M: number to pattern
def number_to_pattern(number, k):
    pattern = ""
    for i in range(k):
        temp = number // 4**(k-i-1)
        if temp == 0:
            pattern += 'A'
        elif temp == 1:
            pattern += 'C'
        elif temp == 2:
            pattern += 'G'
        elif temp == 3:
            pattern += 'T'
        number -= temp * 4**(k-i-1)
    return pattern


# BA01L: pattern_to_number
def pattern_to_number(pattern):
    number = 0
    for i in range(len(pattern)):
        pw = len(pattern) - i - 1
        nuc = 0
        chr = pattern[i]
        if chr == 'A':
            nuc = 0
        elif chr == 'C':
            nuc = 1
        elif chr == 'G':
            nuc = 2
        elif chr == 'T':
            nuc = 3
        pw_val = 4**pw
        number  = number +(pw_val * nuc)
    return number


# BA01L: pattern to number, recursive
def patternToNumber(Pattern):  # solution by hadaarjan (Rosalind site)
    symbolToNumber = {'A':0, 'C':1, 'G':2, 'T':3}
    n = len(Pattern)
    if n == 0:
        return 0
    elif n == 1:
        return symbolToNumber[Pattern]
    else:
        return 4*patternToNumber(Pattern[:-1]) + symbolToNumber[Pattern[-1]]


# BA01M: number to pattern, recursive
def numberToPattern(n, kmer):  # soultion by hadaarjan (Rosalind site)
    numberToSymbol = {0:'A', 1:'C', 2:'G', 3:'T'}
    pattern = ''
    while n > 0:
        remainder = n%4
        pattern = numberToSymbol[remainder] + pattern
        n = n/4
    if kmer - len(pattern) == 0:
        return pattern
    else:
        return (kmer - len(pattern))*'A' + pattern


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01l.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    pattern = lines[0]
    f.close()

    start_time = time.time()
    print(pattern_to_number(pattern))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()