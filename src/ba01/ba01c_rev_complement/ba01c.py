"""
Rosalind: BA1C
Find the Reverse Complement of a String

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C'
and 'G'. Given a nucleotide p, we denote its complementary nucleotide as p. The
reverse complement of a DNA string Pattern = p1…pn is the string rPattern = pn …
p1 formed by taking the complement of each nucleotide in Pattern, then reversing
the resulting string.

For example, the reverse complement of Pattern = "GTCA" is rPattern = "TGAC".

Reverse Complement Problem
Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: rPattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT

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
      -> Find the Most Frequent Words with Mismatches in a String (BA1I)
      -> NEXT: Find Frequent Words with Mismatches and Reverse Complements (BA1)
        -> HERE: Find the Reverse Complement of a String (BA1C)

Plan 1.
- human DNA is double stranded helical structure, and has direction

  5-AAAACCCGGT-3
  3-  ...     -5     I should read this strand by 5 to 3 direction.

- steps:
  a. reverse given DNA strand
  b. get complement nucleotide for each character, that's all
"""
#!/usr/bin/env python
import time


# reverse string
def rev_str(str):
    return str[::-1]


# BA01C: reverse complement
def complement(pattern):
    dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ""
    for nucleotide in pattern:
        complement += dict[nucleotide]
    return complement

# BA01C: reverse complement
def rev_complement(pattern):
    dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ""
    for nucleotide in pattern:
        complement += dict[nucleotide]
    return complement[::-1]

# main function
def main():
    # read lines
    f = open('/home/wsl/rosalind/data/ba01c.txt','r')
    lines = f.readlines()
    pattern = lines[0]
    f.close()

    start_time = time.time()
    print(complement(rev_str(pattern)))
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    print(rev_complement(pattern))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()
