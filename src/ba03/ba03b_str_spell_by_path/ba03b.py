"""
Rosalind: BA3B (difficulty: 1/5)
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern1, ... , Patternn
such that the last k - 1 symbols of Patterni are equal
to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> PREV: Generate the k-mer Composition of a String (BA3A)
      -> HERE: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * NEXT: Construct OverapGraph (BA3C)

TAATGCCATGGGATGTT
TAA->AAT->ATG->TGC->GCC->CCA->CAT->ATG->TGG->GGG->GGA->GAT->ATG->TGT->GTT
"""

import time


def reconstruct(genome_path):
    result = genome_path[0]
    for node in genome_path[1:]:
        result += node[-1]
    return result


def main():
    f = open('/home/wsl/rosalind/data/ba03b.txt', 'r')
    genome_path = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print(reconstruct(genome_path))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()