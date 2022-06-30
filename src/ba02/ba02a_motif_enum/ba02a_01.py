"""
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

    MOTIFENUMERATION(Dna, k, d)
        Patterns <- an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern' differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns

Implanted Motif Problem
Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.

Given: Integers k and d, followed by a collection of strings Dna.

Return: All (k, d)-motifs in Dna.

Sample Dataset
3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output
ATA ATT GTT TTT

═════════════════════════════════════════════════

- by Donny
"""

import time


def hdist(x,y):
    nmm = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            nmm += 1
    return nmm


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


def motifEnumeration (dna, k, d):
    ''' Function to generate all kmers that are present in all strands of dna
        with at most d mismatches '''
    ## Create an empty set for each DNA strand and fill it with all k-mers
    ## with at most d mismatches
    strand_kmers = [set() for _ in range(len(dna))]
    for count,strand in enumerate(dna):
        for i in range(len(strand) - k + 1):
            kmer = strand[i:i+k]
            strand_kmers[count].update(neighbors(kmer,d))
    ## Now intersect all sets to end up with only the sequences that occur
    ## in all strands
    result = strand_kmers[0]
    for i in range(1,len(strand_kmers)):
        result.intersection_update(strand_kmers[i])
    return ' '.join(list(result))


def main():
    f = open('/home/wsl/rosalind/data/ba02a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, d = [int(number.strip()) for number in lines[0].split()]
    strings = [line.strip() for line in lines[1:]]
    f.close()

    start_time = time.time()
    answer = motifEnumeration(strings, k, d)
    print(answer)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()