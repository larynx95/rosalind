"""
Rosalind: BA3A (difficulty: 1/5)
Generate the k-mer Composition of a String

Given a string Text, its k-mer composition Compositionk(Text)
is the collection of all k-mer substrings of Text
(including repeated k-mers).
For example,

Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}

Note that we have listed k-mers in lexicographic order
(i.e., how they would appear in a dictionary)
rather than in the order of their appearance in TATGGGGTGC.
We have done this because the correct ordering of the reads
is unknown when they are generated.

String Composition Problem
Generate the k-mer composition of a string.

Given: An integer k and a string Text.

Return: Compositionk(Text)
        (the k-mers can be provided in any order).

Sample Dataset
5
CAATCCAAC

Sample Output
AATCC
ATCCA
CAATC
CCAAC
TCCAA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> HERE: Generate the k-mer Composition of a String (BA3A)
      -> NEXT: String Reconstruction Problem with k-mer composition(BA3B)

═════════════════════════════════════════════════

References:
- Python Generators <--- use generator if possible!
  https://www.programiz.com/python-programming/generator
- How to loop through a generator
  https://stackoverflow.com/questions/11539194/how-to-loop-through-a-generator
"""

import time


def all_kmers(k, text):
    """ list comprehension version """
    return [text[i:i+k] for i in range(len(text)-k+1)]


def all_kmers_gen(k, text):
    """ generator version """
    for i in range(len(text)-k+1):
        yield text[i:i+k]


def all_kmers_sorted(k, text):
    return sorted(all_kmers(k, text))


def main():
    f = open('/home/wsl/rosalind/data/ba03a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    text = lines[1]
    f.close()

    start_time = time.time()
    kmers = all_kmers(k, text)
    print('\n'.join(kmers))
    print("\t--- %s seconds ---" %(time.time() - start_time))

    start_time = time.time()
    kmers = all_kmers_gen(k, text)
    print('\n'.join(kmers))
    print("\t--- %s seconds ---" %(time.time() - start_time))


if __name__ == "__main__":
    main()