"""
Rosalind: BA1M
Implement n2p

Convert an integer to its corresponding DNA string.

Given: Integers index and k.
Return: n2p(index, k).

Sample Dataset
45
4

Sample Output
AGTC

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
  ║                                                                    ║
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
  ║                                                          ║
  ╚══════════════════════════════════════════════════════════╝
"""

import time


# BA01M: number to pattern
def n2p_ver2(number, k):
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
def p2n_ver2(pattern):
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


# BA01L: pattern to number
def p2n(Pattern):  # solution by hadaarjan (Rosalind site)
    symbolToNumber = {'A':0, 'C':1, 'G':2, 'T':3}
    n = len(Pattern)
    if n == 0:
        return 0
    elif n == 1:
        return symbolToNumber[Pattern]
    else:
        return 4*p2n(Pattern[:-1]) + symbolToNumber[Pattern[-1]]


# BA01M: pattern to number
def n2p(n, kmer):  # soultion by hadaarjan (Rosalind site)
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


# BA01K: compute frequency
def compute_freq(text, k):
    freq_arr = [0] * (4**k)
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        num = p2n_ver2(kmer)
        freq_arr[num] += 1
    return freq_arr


# main
def main():
    f = open('/home/wsl/rosalind/data/ba01m.txt', 'r')
    lines = [line.strip()  for line in f.readlines()]
    number = int(lines[0])
    k = int(lines[1])
    f.close()

    start_time = time.time()
    print(n2p_ver2(number, k))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()