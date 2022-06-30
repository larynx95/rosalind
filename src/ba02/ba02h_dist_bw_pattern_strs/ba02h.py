"""
Rosalind: BA2H
Implement DistanceBetweenPatternAndStrings

The first potential issue with implementing MedianString from "Find a Median String" is
writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai),
the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}.
This task is achieved by the following pseudocode.

  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

Sample Dataset
AAA
TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

Sample Output
5

═════════════════════════════════════════════════

Pseudocode:
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ MEDIANSTRING(Dna, k)                                                   ║
  ║     distance 1                                                         ║
  ║     for i 0 to 4^k - 1                                                 ║
  ║         pattern <- NumberToPattern(i, k)                               ║
  ║         if distance > DistanceBetweenPaternAndStrings(Pattern, Dna)    ║
  ║             distance <- DistanceBetweenPatternAndStrings(Patern, Dna)  ║
  ║             Median <- Pattern                                          ║
  ║     return Median                                                      ║
  ╚════════════════════════════════════════════════════════════════════════╝

"""

from cmath import inf
import time


def all_kmers(dna, k):
    """ get all k-mers from a string """
    return (dna[i:i+k] for i in range(len(dna) - k + 1))


def hdist(str1, str2):
    return sum([str1[i] != str2[i] for i in range(len(str1.strip()))])


def dist_bw_pat_strs(pattern, dnas):
    dist = 0
    k = len(pattern)
    for str in dnas:
        hd = inf
        for kmer in all_kmers(str, k):
            hd_kmer = hdist(pattern, kmer)
            if hd > hd_kmer:
                hd = hd_kmer
        dist += hd
    return dist


def n2p(number, k):
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


def median_string(dnas, k):
    distance = inf
    median = ""
    for i in range(0, 4**k):
        pattern = n2p(i, k)
        if distance > dist_bw_pat_strs(pattern, dnas):
            distance = dist_bw_pat_strs(pattern, dnas)
            median = pattern
    return median


def main():
    f = open('/home/wsl/rosalind/data/ba02h.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    pattern = lines[0]
    dnas = [str.strip() for str in lines[1].split()]
    f.close()

    start_time = time.time()
    print(dist_bw_pat_strs(pattern, dnas))  # --- 0.00015878677368164062 seconds ---
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    print(median_string(dnas, 8))  # --- 2.6731925010681152 seconds ---
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()