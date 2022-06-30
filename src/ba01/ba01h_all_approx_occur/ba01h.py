"""
Rosalind: BA1H
Find All Approximate Occurrences of a Pattern in a String

We say that a k-mer Pattern appears as a substring of Text with at most d
mismatches if there is some k-mer substring Pattern' of Text having d or fewer
mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our
observation that a DnaA box may appear with slight variations leads to the
following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return: All starting positions where Pattern appears as a substring of Text with
at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3

Sample Output
6 7 26 27 78

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
      -> PREV: hamming distance (BA1G)
      -> HERE: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> frequency array
          -> NEXT: number to pattern (BA1M), pattern to number (BA1L)

Plan 1.
- pseudocode:

    APPROXIMATEPATTERNCOUNT(Text, Pattern, d)
      count 0
      for i <- 0 to |Text| - |Pattern|
        Pattern' <- Text(i, |Pattern|)
        if HAMMINGDISTANCE(Pattern, Pattern’) >= d
          count count + 1
      return count

"""

import time


# BA01G: hamming distance
def hdistance(str1, str2):
    count = 0
    for i in range(len(str1.strip())):  # remove the last '\n'
        if str1[i] != str2[i]:
            count += 1
    return count


# BA01H: all approximate occurrences
def approx_occur(pattern, genome, d):
    indices = []
    for i in range(len(genome) - len(pattern) + 1):
        subs = genome[i:i+len(pattern)]
        if hdistance(pattern, subs) <= d:
            indices.append(i)
    return indices


# BA01H: all approximate occurrences
def approx_pattern_count(pattern, genome, d):
    count = 0
    for i in range(len(genome) - len(pattern) + 1):
        subs = genome[i:i+len(pattern)]
        if hdistance(pattern, subs) <= d:
            count += 1
    return count


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01h.txt','r')
    lines = f.readlines()
    pattern = lines[0].strip()
    genome = lines[1].strip()
    d = int(lines[2].strip())
    f.close()

    start_time = time.time()
    print(approx_occur(pattern, genome, d))
    print(len(approx_occur(pattern, genome, d)))
    print(approx_pattern_count(pattern, genome, d))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()