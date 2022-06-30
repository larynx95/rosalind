"""
Rosalind: BA1K
Generate the Frequency Array of a String

Given an integer k, we define the frequency array of a string Text as an array
of length 4k, where the i-th element of the array holds the number of times that
the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

kmer      AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
index      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
frequency  3   0   2   0   1   0   0   0   0   1   3   1   0   0   1   0

Computing a Frequency Array
Generate the frequency array of a DNA string.

Given: A DNA string Text and an integer k.

Return: The frequency array of k-mers in Text.

Sample Dataset
ACGCGGCTCTGAAA
2

Sample Output
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

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
      -> Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - by frequency array
        -> HERE: frequency array
          -> PREV: number to pattern (BA1M), pattern to number (BA1L)
        -> NEXT: neigbors (BA1N)

Plan 1.
- pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  COMPUTINGFREQUENCIES(Text, k)                    ║
  ║    for i <- 0 to 4^k- 1                           ║
  ║      FREQUENCYARRAY(i) <- 0                       ║
  ║    for i <- 0 to |Text| - k                       ║
  ║      Pattern <- Text(i, k)                        ║
  ║      j <- PATTERNTONUMBER(Pattern)                ║
  ║      FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1   ║
  ║    return FREQUENCYARRAY                          ║
  ╚═══════════════════════════════════════════════════╝

"""

import time


# BA01L: pattern to number
def pattern_to_number(pattern):
    """
    str -> int
    returns an integer from a pattern
    """
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


# BA01K: compute frequency
def compute_freq(text, k):
    """
    (str,int) -> [int]
    >>> compute_freq('ACGCGGCTCTGAAA', 2)
        [2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]
    """
    freq_arr = [0] * (4**k)
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        num = pattern_to_number(kmer)
        freq_arr[num] += 1
    return freq_arr


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01k.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    pattern = lines[0]
    k = int(lines[1])
    f.close()

    # check time
    start_time = time.time()
    print(compute_freq(pattern, k))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()