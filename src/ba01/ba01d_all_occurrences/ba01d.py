"""
Rosalind: BA1D
Find All Occurrences of a Pattern in a String

In this problem, we ask a simple question: how many times can one string occur
as a substring of another? Recall from “Find the Most Frequent Words in a
String” that different occurrences of a substring can overlap with each other.
For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
Find all occurrences of a pattern in a string.

Given: Strings Pattern and Genome.

Return: All starting positions in Genome where Pattern appears as a substring.
Use 0-based indexing.

Sample Dataset
ATAT
GATATATGCATATACTT

Sample Output
1 3 9

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
      -> HERE: find all occurrence of a pattern in a string (BA1D)
      -> NEXT: Clump finidng problem (BA1E)
"""
#!/usr/bin/env python
import time


# BA01D: find all occurrences
def find_all_occur(pattern, dna):
    """
    (str,str) -> [int]
    >>> all_occurrence('ATAT','GATATATGCATATACTT')
        [1,3,9]
    """
    indices = []
    for i in range(0, len(dna)-len(pattern)+1):
        if pattern == dna[i:i+len(pattern)]:
            indices.append(i)
    return indices


# main functino
def main():
    f = open('/home/wsl/rosalind/data/ba01d.txt', 'r')
    lines = f.readlines()
    pattern = lines[0].strip()
    dna = lines[1].strip()
    f.close()

    start_time = time.time()
    print(find_all_occur(pattern, dna))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()