"""
Rosalind: BA1J
Find Frequent Words with Mismatches and Reverse Complements

We now extend "Find the Most Frequent Words with Mismatches in a String" to find
frequent words with both mismatches and reverse complements.
Recall that rPattern refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
Find the most frequent k-mers (with mismatches and reverse complements) in a DNA
string.

Given: A DNA string Text as well as integers k and d.

Return: All k-mers Pattern maximizing the sum
        Count_d(Text, Pattern) + Count_d(Text, rPattern)
        over all possible k-mers.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
ATGT ACAT

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
      -> HERE: Find Frequent Words with Mismatches and Reverse Complements (BA1J)
        -> Find the Reverse Complement of a String (BA1C)
      ↓
    * Solving the most frequent words with mismatch by improving algorithms
      -> Generate the Frequency Array of a String (BA1K)
        -> NEXT: pattern to number (BA1L), number to pattern (BA1M)

"""

import time


# BA01C: reverse complement
def rev_complement(pattern):
    dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ""
    for nucleotide in pattern:
        complement += dict[nucleotide]
    return complement[::-1]


# BA01G: hamming distance
def hdist(str1, str2):
    count = 0
    for i in range(len(str1.strip())):  # remove the last '\n'
        if str1[i] != str2[i]:
            count += 1
    return count


# BA01L: pattern to number
def pattern_to_number(Pattern):  # soultion by hadaarjan (Rosalind site)
    symbolToNumber = {'A':0, 'C':1, 'G':2, 'T':3}
    n = len(Pattern)
    if n == 0:
        return 0
    elif n == 1:
        return symbolToNumber[Pattern]
    else:
        return 4*pattern_to_number(Pattern[:-1]) + symbolToNumber[Pattern[-1]]


# BA01M: number to pattern
def number_to_pattern(n, kmer):  # soultion by hadaarjan (Rosalind site)
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
def freq_words_mismatch_rc(text, k, d):
    """
    (str,int,int) -> {str}
    returns most frequent word with mismatches reverse complement
    """
    words = set()
    count = {}
    max_count = 0
    for i in range(len(text)-k+1):
        nbs = neighbors(text[i:i+k], d)
        for nb in nbs:
            idx = pattern_to_number(nb)
            if idx not in count:
                count[idx] = 1
            else:
                count[idx] += 1
                if count[idx] > max_count:
                    max_count = count[idx]
        nbs_rc = neighbors(rev_complement(text[i:i+k]), d)
        for nb in nbs_rc:
            idx = pattern_to_number(nb)
            if idx not in count:
                count[idx] = 1
            else:
                count[idx] += 1
                if count[idx] > max_count:
                    max_count = count[idx]
    for cnt in count:
        if count[cnt] == max_count:
            words.add(number_to_pattern(cnt, k))
    return words


# main
def main():
    f = open('/home/wsl/rosalind/data/ba01j.txt', 'r')
    lines = f.readlines()
    text = lines[0].strip()
    k, d = [int(num.strip()) for num in lines[1].split()]
    f.close()


    start_time = time.time()
    answers = freq_words_mismatch_rc(text, k, d)
    for answer in answers:
        print(answer, end=" ")
    print()
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()

