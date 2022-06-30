"""
Rosalind: BA1I
Find the Most Frequent Words with Mismatches in a String

We defined a mismatch in "Compute the Hamming Distance Between Two Strings". We
now generalize "Find the Most Frequent Words in a String" to incorporate
mismatches as well.

Given strings Text and Pattern as well as an integer d, we define Countd(Text,
Pattern) as the total number of occurrences of Pattern in Text with at most d
mismatches. For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because
AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA,
AAACA, and AAAGA. Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in Text is simply a string Pattern
maximizing Count_d(Text, Pattern) among all k-mers. Note that Pattern does not
need to actually appear as a substring of Text; for example, AAAAA is the most
frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA
does not appear exactly in this string. Keep this in mind while solving the
following problem.

Frequent Words with Mismatches Problem
Find the most frequent k-mers with mismatches in a string.

Given: A string Text as well as integers k and d.

Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
GATG ATGC ATGT

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
      -> HERE: Find the Most Frequent Words with Mismatches in a String (BA1I)
        -> frequency array
          -> number to pattern (BA1M), pattern to number (BA1L)
        -> PREV: neigbors (BA1N)

Plan 1.
- Pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  APPROXIMATEPATTERNCOUNT(Text, Pattern, d)        ║
  ║    count 0                                        ║
  ║    for i 0 to |Text| - |Pattern|                  ║
  ║      Pattern <- Text(i, |Pattern|)                ║
  ║      if HAMMINGDISTANCE(Pattern, Pattern’) >= d   ║
  ║        count count + 1                            ║
  ║    return count                                   ║
  ╚═══════════════════════════════════════════════════╝


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


  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  FREQUENTWORDSWITHMISMATCHES(Text, k, d)                              ║
  ║    FrequentPatterns an empty set                                      ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLOSE(i) <- 0                                                    ║
  ║      FREQUENCYARRAY <- 0                                              ║
  ║    for i <- 0 to |Text| - k                                           ║
  ║      Neighborhood <- NEIGHBORS(Text(i, k), d)                         ║
  ║      for each Pattern from Neighborhood                               ║
  ║      index <- p2n(Pattern)                                            ║
  ║      CLOSE(index) <- 1                                                ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLOSE(i) = 1                                                  ║
  ║      Pattern <- n2p(i, k)                                             ║
  ║      FREQUENCYARRAY(i) <- APPROXIMATEPATTERNCOUNT(Text, Pattern, d)   ║
  ║    maxCount maximal value in FREQUENCYARRAY                           ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if FREQUENCYARRAY(i) = maxCount                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝


  ╔══════════════════════════════════════════════════════════════════════════╗
  ║  FINDINGFREQUENTWORDSWITHMISMATCHESBYSORTING(Text, k, d)                 ║
  ║    FrequentPatterns <- an empty set                                      ║
  ║    Neighborhoods <- an empty list                                        ║
  ║    for i <- 0 to |Text| - k                                              ║
  ║      add NEIGHBORS(Text(i, k), d) to Neighborhoods                       ║
  ║    form an array NEIGHBORHOODARRAY holding all strings in Neighborhoods  ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      Pattern <- NEIGHBORHOODARRAY(i)                                     ║
  ║      INDEX(i) <-  p2n(Pattern)                                           ║
  ║      COUNT(i) <- 1                                                       ║
  ║    SORTEDINDEX SORT(INDEX)                                               ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      if SORTEDINDEX(i) = SORTEDINDEX(i + 1)                              ║
  ║        COUNT(i + 1) <- COUNT(i) + 1                                      ║
  ║    maxCount maximum value in array COUNT                                 ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║     if COUNT(i) = maxCount                                               ║
  ║       Pattern <- n2p(SORTEDINDEX(i), k)                                  ║
  ║       add Pattern to FrequentPatterns                                    ║
  ║    return FrequentPatterns                                               ║
  ╚══════════════════════════════════════════════════════════════════════════╝

"""

import time


#################################################
# brute-force algorithm (inefficient)
#################################################


def all_kmers(characters, k):
    import itertools as it
    return [''.join(c) for c in it.product(characters, repeat=k)]


def all_kmers_iter(k, y=''):
    if k == 0:
        yield y
    else:
        for m in ['A', 'C', 'T', 'G']:
            yield from all_kmers_iter(k - 1, m + y)


# BA01I: frequent words mismatches
def freq_words_mismatch_brute(dna, k, d):
    allKmers = all_kmers_iter(k)
    tempDic = dict()
    answer = []
    maxCount = 0
    for s in allKmers:
        for i in range(len(dna) - k + 1):
            if hdist(s, dna[i:i + k]) <= d:
                if s in tempDic:
                    tempDic[s] += 1
                else:
                    tempDic[s] = 1
    for k in tempDic.keys():
        if tempDic[k] > maxCount:
            maxCount = tempDic[k]
    for k in tempDic.keys():
        if tempDic[k] == maxCount:
            answer.append(k)
    return answer


#################################################
# using frequency array
#################################################

# BA01M: number to pattern
def num_to_pattern(number, k):
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


# BA01L: pattern to number
def pattern_to_num(pattern):
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
    freq_arr = [0] * (4**k)
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        num = pattern_to_num(kmer)
        freq_arr[num] += 1
    return freq_arr


# BA01G: hamming distance
def hdist(str1, str2):
    count = 0
    for i in range(len(str1.strip())):  # remove the last '\n'
        if str1[i] != str2[i]:
            count += 1
    return count


# BA01G: hamming distance
def hdist_lst_comprehension(x, y):
    return sum([x[i] != y[i] for i in range(len(x))])


# BA01H: all approximate occurrences
def approx_pattern_count(pattern, genome, d):
    count = 0
    for i in range(len(genome) - len(pattern) + 1):
        subs = genome[i:i+len(pattern)]
        if hdist(pattern, subs) <= d:
            count += 1
    return count


# BA01L: pattern to number
def p2n(Pattern):  # soultion by hadaarjan (Rosalind site)
    symbolToNumber = {'A':0, 'C':1, 'G':2, 'T':3}
    n = len(Pattern)
    if n == 0:
        return 0
    elif n == 1:
        return symbolToNumber[Pattern]
    else:
        return 4*p2n(Pattern[:-1]) + symbolToNumber[Pattern[-1]]


# BA01M: number to pattern
def n2p(n, kmer):  # soultion by hadaarjan (Rosalind site)
    numberToSymbol = {0:'A', 1:'C', 2:'G', 3:'T'}
    pattern = ''
    while n > 0:
        remainder = n % 4
        pattern = numberToSymbol[remainder] + pattern
        n = n//4
    if kmer - len(pattern) == 0:
        return pattern
    else:
        return (kmer - len(pattern))*'A' + pattern


# BA01N: immediate neighbors
def immediate_neighbors(pattern):
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


# BA01N: immediate neighbors
def immediate_neighbors_set(pattern):
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


# BA01N: immdediate neighbors
def immediate_neighbors_set2(pattern):
    results = {pattern}
    for i in range(len(pattern)):
        str_rest = pattern[:i] + pattern[i+1:]
        ch_focused = pattern[i]
        for nuc in "ACGT":
            if ch_focused != nuc:
                results.add((str_rest[:i] + nuc + str_rest[i:]))
    return results


# BA01N: neighbors
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


# BA01I: frequent words mismatches
def freq_words_mismatch(text, k, d):
    """
    (str,int,int) -> {str}
    soultion by hadaarjan (Rosalind site) - algorithm in textbook
    """
    frequentPatterns = set()
    count = [0]*((4**k))
    for i in range(len(text) -k + 1):
        neighborhood = neighbors(text[i:i+k], d)
        for word in neighborhood:
            index = p2n(word)
            count[index] += 1
    maxCount = max(count)
    for i in range(4**k):
        if count[i] == maxCount:
            pattern = n2p(i,k)
            frequentPatterns.add(pattern)
    return frequentPatterns


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01i.txt', 'r')
    lines = f.readlines()
    text = lines[0].strip()
    k, d = [int(num.strip()) for num in lines[1].split()]
    f.close()

    start_time = time.time()
    print(*freq_words_mismatch_brute(text, k, d))
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    print(*freq_words_mismatch(text, k, d))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()
