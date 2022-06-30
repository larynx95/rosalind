"""
Rosalind: BA1F
Find a Position in a Genome Minimizing the Skew

Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference
between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi
(Genome) denote the prefix (i.e., initial substring) of Genome of length i. For
example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i
(from 0 to |Genome|).

Sample Dataset
CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

Sample Output
53 97

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> HERE: minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * NEXT: Most frequent words problem

"""
#!/usr/bin/env python
import time


# BA01F: min skew
def min_skew(dna):
    # get skew array
    arr = [0]
    skew = 0
    for nuc in dna:
        if nuc == 'C':
            skew -= 1
        if nuc == 'G':
            skew += 1
        arr.append(skew)

    # find minimum skew
    min_skew = 0
    for i in range(len(arr)):
        if arr[i] < min_skew:
            min_skew = arr[i]

    # get indices of element with minimum skewness
    indices = []
    for i in range(len(arr)):
        if arr[i] == min_skew:
            indices.append(i)
    return indices


# main
def main():
    f = open('/home/wsl/rosalind/data/ba01f.txt', 'r')
    lines = f.readlines()
    genome = lines[0].strip()
    f.close()

    start_time = time.time()
    print(min_skew(genome))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()