"""
Rosalind: BA6E
Find All Shared k-mers of a Pair of Strings

We say that a k-mer is shared by two genomes
if either the k-mer or its "reverse complement" appears in each genome.
In Figure 1 are four pairs of 3-mers that are shared by "AAACTCATC" and "TTTCAAATC".

A shared k-mer can be represented by an ordered pair (x, y),
where x is the starting position of the k-mer in the first genome
and y is the starting position of the k-mer in the second genome.
For the genomes "AAACTCATC" and "TTTCAAATC",
these shared k-mers are (0,4), (0,0), (4,2), and (6,6).

Shared k-mers Problem
Given two strings, find all their shared k-mers.

Given:
An integer k and two strings.

Return:
All k-mers shared by these strings, in the form of ordered pairs (x, y)
corresponding to starting positions of these k-mers in the respective strings.

Sample Dataset
3
AAACTCATC
TTTCAAATC

Sample Output
(0, 4)
(0, 0)
(4, 2)
(6, 6)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> Cycle to Chromosome (BA6G)
      -> Genome to ColoredEdges (BA6H)
      -> ColoredEdges to Genome (BA6I)
      -> 2-Break on ColoredEdges (BA6J)
      -> 2-Break on Genome (BA6K)
      ↓
    * 2-Break distance (BA6C)
      ↓
    * PREV: 2-Break sorting (BA6D)
      ↓
    * HERE: Shared k-mers (BA6E)

Info.
  * dot plot (fig 6.20)
    AGCAGGTTATCTCCCTGT, K=2                  AGCAGGTTATCTCCCTGT, K=3
    T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    G ─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─  G ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─  T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─
    C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─O─┼─┼─┼─  C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─
    C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─  C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─
    C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─  C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─
    T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─O─┼─┼─┼─┼─┼─┼─  T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─
    C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─O─┼─┼─┼─  C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─
    T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─O─┼─┼─┼─┼─┼─┼─  T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─
    A ─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─  A ─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T ─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  T ─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T ─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  T ─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    G ─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─  G ─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    G ─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  G ─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A ─O─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  A ─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    C ─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  C ─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    G ─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  G ─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A ─O─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  A ─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
       A G C A G G T T A T C T C C C T G T      A G C A G G T T A T C T C C C T G T

    AGCAGGTTATCTCCCTGT, K=3, rev complement    AGCAGGTTATCTACCTGT & AGCAGGAGATAAACCTGT, k=3, rev complement
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    C-G ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  C-G ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─  A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─
    G-C ─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─  G-C ─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─
    G-C ─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─  G-C ─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─
    G-C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─  T-A ─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─ /T-A ─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─ reversal TTATCT
    G-C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─ /T-A ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─ /A-T ─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T-A ─┼─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─ /T-A ─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─ /C-G ─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─ /T-A ─┼─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─┼─┼─┼─┼─
    C-G ─┼─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  C-G ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    C-G ─┼─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  C-G ─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T-A ─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─  T-A ─┼─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─┼─
    G-C ─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─  G-C ─┼─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─R─┼─┼─┼─
    C-G ─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  C-G ─┼─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T-A ─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─  T-A ─O─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─ TODO: What does this mean?
         A G C A G G[T T A T C T C C C T G T        A G C A G G/T T A T C T/A C C T G T        Is this just for visualization?
         | | | | | | | | | | | | | | | | | |        | | | | | | | | | | | | | | | | | |
         T C G T C C A A T A G A G G G A C A        T C G T C C/A A T A G A/T G G A C A

  * sample dataset
           k=3
    idx    012          012              456            678
    str1   AAACTCATC    AAACTCATC    AAACTCATC    AAACTCATC
    rev    TTT          TTT              AGT            TAG
    str2   TTTCAAATC    TTTCAAATC    TTTCAAATC    TTTCAAATC
    rev        TTT      AAA            AGT              TAG
    idx        456      012            234              678

    G-C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T-A ─┼─┼─┼─┼─┼─┼─O─┼─┼─
    T-A ─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    T-A ─O─┼─┼─┼─┼─┼─┼─┼─┼─
    G-C ─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─O─┼─┼─┼─┼─
    A-T ─┼─┼─┼─┼─┼─┼─┼─┼─┼─
    A-T ─R─┼─┼─┼─┼─┼─┼─┼─┼─
         A A A C T C A T C
         | | | | | | | | |
         T T T G A G T A G

Plan 1. (Failed - too slow for long strings)
  * simple for-loops
  * for k-mer in string1, find the same k-mer or complement k-mer in string2

Plan 2.
  * modifying frequency array in chapter 01
    number to pattern (BA1M), pattern to number (BA1L), frequency array (BA1K)
  * plus, using faster data struncture (e.g. dictionary, not list)

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time
import sys


def complement(pattern):
    """
    str -> str
    returns a complement
    """
    dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ""
    for nucleotide in pattern:
        complement += dict[nucleotide]
    return complement


def reverse_complement(pattern):
    """
    str -> str
    returns a reverse complement
    >>> reverse_complement('ACTTG')
        CAAGT
    """
    dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ""
    for nucleotide in pattern:
        complement += dict[nucleotide]
    return complement[::-1]


def shared_kmer_slow(k, astr, bstr):
    """
    (int,str,str) -> [(int,int)]
    returns a list of shared k-mers in the form of ordered pairs (x, y)
    """
    pairs = []
    for i in range(len(astr)-k+1):
        kmer_a = astr[i:i+k]
        for j in range(len(bstr)-k+1):
            kmer_b = bstr[j:j+k]
            if kmer_a == kmer_b or kmer_a == reverse_complement(kmer_b):
                pairs.append((i,j))
    return pairs


def num_to_pattern(number, k):
    """
    (int,int) -> str
    returns a pattern from an integer (BA1M)
    """
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


def pattern_to_num(pattern):
    """
    str -> int
    returns an integer from a pattern (BA1L)
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


def compute_freq(text, k):
    """
    (str,int) -> [int]
    returns a frequency array (BA1K)
    >>> compute_freq('ACGCGGCTCTGAAA', 2)
        [2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]
    """
    freq_arr = [0] * (4**k)
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        num = pattern_to_num(kmer)
        freq_arr[num] += 1
    return freq_arr


def freq_dic(text, k):
    """
    (str,int) -> {int:[int]}
    returns a frequency array as python dictionary
    key: number representing a k-mer pattern
    value: indices of k-mer and its reverse complement
    >>> freq_dic('AAAACTC', 3)
        {0:[0,1],63:[0,1],1:[2],47:[2],7:[3],11:[3],29:[4],34:[4]}
    """
    freq = {}
    for i in range(len(text)-k+1):
        kmer    = text[i:i+k]               # kmer
        num_kmer    = pattern_to_num(kmer)
        if num_kmer not in freq:
            freq[num_kmer] = [i]
        else:
            freq[num_kmer].append(i)
        kmer_rc = reverse_complement(kmer)  # k-mer reverse complement
        num_kmer_rc = pattern_to_num(kmer_rc)
        if num_kmer_rc not in freq:
            freq[num_kmer_rc] = [i]
        else:
            freq[num_kmer_rc].append(i)
    return freq


def shared_kmer(k, astr, bstr):
    """
    (int,str,str) -> [(int,int)]
    returns a list of shared k-mers in the form of ordered pairs (x, y)
    """
    result = []
    freq = freq_dic(bstr, k)
    for i in range(len(astr)-k+1):
        kmer = astr[i:i+k]
        num = pattern_to_num(kmer)
        if num in freq:
            for j in freq[num]:
                result.append((i,j))
    return result


def main():
    f = open('/home/wsl/rosalind/data/ba06e.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    astr = lines[1]
    bstr = lines[2]
    f.close()

    start_time = time.time()
    pairs = shared_kmer(k,astr, bstr)
    # if too many lines to display on terminal
    if len(pairs) > 100:
        with open('/home/wsl/rosalind/data/ba06e_output.txt', 'w') as f:
            original_stdout = sys.stdout
            sys.stdout = f
            for pair in pairs:
                print(pair)
            sys.stdout = original_stdout
    else:
        for pair in pairs:
            print(pair)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
