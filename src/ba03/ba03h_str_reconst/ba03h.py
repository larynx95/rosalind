"""
Rosalind: BA3H
Reconstruct a String from its k-mer Composition

String Reconstruction Problem
Reconstruct a string from its k-mer composition.

Given: An integer k followed by a list of k-mers Patterns.

Return: A string Text with k-mer composition equal to Patterns.
(If multiple answers exist, you may return any one.)

Sample Dataset
4
CTTA
ACCA
TACC
GGCT
GCTT
TTAC

Sample Output
GGCTTACCA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * PREV: Construct OverapGraph (BA3C)
      -> HERE: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> NEXT: Construct De Bruijn Graph (BA3D)

Plan 1.
  * kmers -> Overlap graph -> Hamiltonian path
  * split k-mer to two (k-1)-mers (prefix and suffix)
    CTTA -> CTT,TTA
    ACCA -> ACC,CCA
    TACC -> TAC,ACC
    GGCT -> GGC,GCT
    GCTT -> GCT,CTT
    TTAC -> TTA,TAC
    prefixes: [CCT,ACC,TAC,GGC,GCT,TTA]
    suffixes: [TTA,CCA,ACC,GCT,CTT,TAC]

Plan 2.
  * kmers -> De Bruijn graph -> Eulerian path

═════════════════════════════════════════════════

References:
- Why can’t you modify lists through "for in" loops in Python?
  https://www.quora.com/Why-can%E2%80%99t-you-modify-lists-through-for-in-loops-in-Python
"""
#!/usr/bin/env python
import time


def str_recon_wrong_01(kmers):
    """
    TODO: Fix it! problem in modifying list in for-loop
    """
    pattern = kmers.pop(0)
    for kmer in kmers:
        if pattern[-3:] == kmer[:-1]:
            pattern += kmer[-1]
            kmers.remove(kmer)
        if not kmers:
            return pattern
        if pattern[:3] == kmer[1:]:
            pattern = kmer[0] + pattern
            kmers.remove(kmer)
    return pattern


def str_recon_wrong_02(kmers):
    """
    TODO: Fix it! problem in modifying list in while-loop
    """
    pattern = kmers.pop(0)
    iter_kmer = iter(kmers)
    while kmers:
        kmer = next(iter_kmer)
        if pattern[-3:] == kmer[:-1]:
            pattern += kmer[-1]
            kmers.remove(kmer)
        if pattern[:3] == kmer[1:]:
            pattern = kmer[0] + pattern
            kmers.remove(kmer)
    return pattern


def str_recon(k, kmers):
    pattern = kmers.pop(0)        # starting kmer for string enlongation
    i = 0                         # variable for adjusting index of list
    while kmers:
        kmer = kmers[i]
        i += 1
        # if there's a consecutive kmer -> appending process
        if pattern[-(k-1):] == kmer[:k-1]:
            pattern += kmer[-1]
            kmers.remove(kmer)
            i = 0
        # if there's no consecutive kmer -> prepending process
        if pattern[:k-1] == kmer[-(k-1):]:
            pattern = kmer[0] + pattern
            kmers.remove(kmer)
            i = 0
    return pattern


def str_recon_recur(k, kmers):
    """ TODO: Implement recursive fuction """
    pass


def main():
    f = open('/home/wsl/rosalind/data/ba03h.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    kmers = lines[1:]
    f.close()

    # k_sample = 4
    # kmers_sample = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
    # print(str_recon(k_sample, kmers_sample))

    start_time = time.time()
    print(str_recon(k, kmers))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
