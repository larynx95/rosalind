"""
Rosalind: BA1N
Generate the d-Neighborhood of a String

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming
distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String
Find all the neighbors of a pattern.

Given: A DNA string Pattern and an integer d.
Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
ACG
1

Sample Output
CCG
TCG
GCG
AAG
ATG
AGG
ACA
ACC
ACT
ACG

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
      -> NEXT: Find the Most Frequent Words with Mismatches in a String (BA1I)
        -> PREV: frequency array
          -> number to pattern (BA1M), pattern to number (BA1L)
        -> HERE: neigbors (BA1N)

Plan 1.
- pseudocode
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
  ║    Neighborhood <- the set consisting of the single string Pattern    ║
  ║    for i = 1 to |Pattern|                                             ║
  ║      symbol <- i-th nucleotide of Pattern                             ║
  ║      for each nucleotide x different from symbol                      ║
  ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
  ║        add Neighbor to Neighborhood                                   ║
  ║    return Neighborhood                                                ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NEIGHBORS(Pattern, d)                                   ║
  ║    if d = 0                                              ║
  ║      return {Pattern}                                    ║
  ║    if |Pattern| = 1                                      ║
  ║      return {A, C, G, T}                                 ║
  ║    Neighborhood <- an empty set                          ║
  ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
  ║    for each string Text from SuffixNeighbors             ║
  ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
  ║        for each nucleotide x                             ║
  ║          add x + Text to Neighborhood                    ║
  ║      else                                                ║
  ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
  ║    return Neighborhood                                   ║
  ╚══════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════╗
  ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
  ║    Neighborhood <- set consisting of single string Pattern   ║
  ║    for j = 1 to d                                            ║
  ║      for each string Pattern' in Neighborhood                ║
  ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
  ║        remove duplicates from Neighborhood                   ║
  ║    return Neighborhood                                       ║
  ╚══════════════════════════════════════════════════════════════╝
"""

import time


# BA01G: hamming distance
def hdist(x,y):
    nmm = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            nmm += 1
    return nmm


# BA01N: neighbors, immediate
def immediate_neighbors(pattern):
    """
    str -> [str]
    >>> immediate_neighbors("ACG")
        ['CCG','GCG','TCG','AAG','AGG','ATG','ACA','ACC','ACT','ACG']
    """
    results = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        nucleotides = "ACGT"
        for chr in nucleotides:
            if symbol != chr:
                pre = pattern[0:i]
                post = pattern[i+1: len(pattern)-i+1]
                temp = pre + chr + post
                results.append(temp)
    results.append(pattern)
    return results


# BA01N: neighbors, immediate
def immediate_neighbors_set(pattern):
    """
    str -> {str}
    >>> immediate_neighbors_set("ACG")
        {'CCG','ACA','ACG','AGG','AAG','ACC','GCG','ATG','ACT','TCG'}
    """
    results = set()
    for i in range(len(pattern)):
        symbol = pattern[i]
        nucleotides = "ACGT"
        for chr in nucleotides:
            if symbol != chr:
                pre = pattern[0:i]
                post = pattern[i+1: len(pattern)-i+1]
                temp = pre + chr + post
                results.add(temp)
    results.add(pattern)
    return results


# BA01N: neighbors, immediate
def immediate_neighbors_set2(pattern):
    """
    str -> {str}
    >>> immediate_neighbors_set2("ACG")
        {'AGG','ACT','CCG','ACC','ACG','GCG','ATG','AAG','ACA','TCG'}
    """
    results = {pattern}
    for i in range(len(pattern)):
        str_rest = pattern[:i] + pattern[i+1:]
        ch_focused = pattern[i]
        for nuc in "ACGT":
            if ch_focused != nuc:
                results.add((str_rest[:i] + nuc + str_rest[i:]))
    return results


# BA01N: neighbors, iterative
def neighbors_iter(pattern, d):
    """
    (str,int) -> {str}
    iterative neighbors fx
    >>> neighbors_iter('ACG',1)
        {'TCG','CCG','GCG','ACG','AAG','AGG','ACC','ATG','ACT','ACA'}
    """
    neighborhood = {pattern}
    for j in range(d):
        for pat in neighborhood.copy():  # <-- set can't be changed during iteration
            im_neighbors = immediate_neighbors_set(pat)
            neighborhood.update(im_neighbors)
    return neighborhood


# BA01N: neighbors, recursive
def neighbors(pattern, d):
    """
    (str,int) -> {str}
    returns all k-mers of Hamming distance at most d from Pattern.
    >>> neighbors('ACG',1)  # <-- this is immediate neighbors
        {'ACC','CCG','AGG','AAG','GCG','ATG','ACT','ACA','TCG','ACG'}
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 0:
        return {}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffixneighbors = neighbors(pattern[1:], d)
    for text in suffixneighbors:
        if hdist(text, pattern[1:]) < d:
            for nt in ['A', 'C', 'G', 'T']:
                neighborhood.add(nt + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01n.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    pattern = lines[0]
    d = int(lines[1])
    f.close()

    start_time = time.time()
    for n in neighbors(pattern, d):
        print(n)
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()
