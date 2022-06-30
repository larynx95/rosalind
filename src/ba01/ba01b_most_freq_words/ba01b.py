"""
Rosalind: BA1B
Find the Most Frequent Words in a String

We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text,
Pattern) among all k-mers. For example, "ACTAT" is a most frequent 5-mer in
"ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most frequent 3-mer of
"CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4

Sample Output
CATG GCAT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> HERE: find frequent words in a string (BA1B)
      -> NEXT: find all occurrence of a pattern in a string (BA1D)

Plan 1.
- writing down
  kmer    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
  index    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  x % 4    0  1  2  3  0  1  2  3  0  1  2  3  0  1  2  3
  x / 4    0  0  0  0  1  1  1  1  2  2  2  2  3  3  3  3

                A:0, C:1, G:2, T:3

    exp   1 0              2 1 0
    Nuc   T T : 3 3        G A T : 2 0 3
          │ └── 4^0*3      │ │ └── 4^0*2
          └──── 4^1*3      │ └──── 4^1*0
                 = 15      └────── 4^2*3
                                    = 50
- Pseudocode:
╔═════════════════════════════════════════════════╗
║ COMPUTINGFREQUENCIES(Text, k)                   ║
║   for i <- 0 to 4^k- 1                          ║
║     FREQUENCYARRAY(i) <- 0                      ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- Text(i, k)                       ║
║     j <- PATTERNTONUMBER(Pattern)               ║
║     FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1  ║
║   return FREQUENCYARRAY                         ║
╚═════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════╗
║ FREQUENTWORDS(Text, k)                          ║
║   FrequentPatterns an empty set                 ║
║     for i <- 0 to |Text| - k                    ║
║     Pattern <- the k-mer Text(i, k)             ║
║     COUNT(i) <- PATTERNCOUNT(Text, Pattern)     ║
║   maxCount maximum value in array COUNT         ║
║   for i <- 0 to |Text| - k                      ║
║     if COUNT(i) = maxCount                      ║
║       add Text(i, k) to FrequentPatterns        ║
║   remove duplicates from FrequentPatterns       ║
║   return FrequentPatterns                       ║
╚═════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════╗
║ FASTERFREQUENTWORDS(Text , k)                     ║
║   FrequentPatterns <- an empty set                ║
║   FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k) ║
║   maxCount <- maximal value in FREQUENCYARRAY     ║
║   for i <- 0 to 4k - 1                            ║
║     if FREQUENCYARRAY(i) = maxCount               ║
║       Pattern <- NUMBERTOPATTERN(i, k)            ║
║       add Pattern to the set FrequentPatterns     ║
║   return FrequentPatterns                         ║
╚═══════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════╗
║ FINDINGFREQUENTWORDSBYSORTING(Text , k)              ║
║   FrequentPatterns <- an empty set                   ║
║   for i <- 0 to |Text| - k                           ║
║     Pattern <- Text(i, k)                            ║
║     INDEX(i) <- PATTERNTONUMBER(Pattern)             ║
║     COUNT(i) <- 1                                    ║
║   SORTEDINDEX <- SORT(INDEX)                         ║
║   for i <- 1 to |Text| - k                           ║
║     if SORTEDINDEX(i) = SORTEDINDEX(i - 1)           ║
║       COUNT(i) = COUNT(i - 1) + 1                    ║
║   maxCount <- maximum value in the array COUNT       ║
║   for i <- 0 to |Text| - k                           ║
║     if COUNT(i) = maxCount                           ║
║       Pattern <- NUMBERTOPATTERN(SORTEDINDEX(i), k)  ║
║       add Pattern to the set FrequentPatterns        ║
║   return FrequentPatterns                            ║
╚══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- Create an empty list in Python with certain size
  https://stackoverflow.com/questions/10712002/create-an-empty-list-in-python-with-certain-size
- Using Python's max to return two equally large values
  https://stackoverflow.com/questions/9853302/using-pythons-max-to-return-two-equally-large-values
"""
#!/usr/bin/env python
import time


# BA02: frequent words
def freq_words(dna, k):
    """
    (str,int) -> {str}
    >>> freq_words('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
        {'GCAT', 'CATG'}
    """
     # construct a dictionary, get the max value
    dict = {}
    max_val = 0
    for i in range(0, len(dna)-k+1):
        kmer = dna[i:i+k]
        if kmer not in dict:
           dict[kmer] = 1
        else:
            dict[kmer] += 1
            if dict[kmer] > max_val:
               max_val = dict[kmer]
    # get a result set of all keys with max value
    result = set()
    for k, v in dict.items():
        if v == max_val:
            result.add(k)
    return result


# BA01M: number to pattern
def number_to_pattern(number, k):
    """
    (int,int) -> str
    >>> number_to_pattern(10, 3)
        AGG
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


# BA01L: pattern to number
def pattern_to_number(pattern):
    """
    (str) -> int
    >>> pattern_to_number('ACGT')
        27
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


# BA01K: compute frequency
def compute_freq(text, k):
    """
    (str,int) -> [int]
    >>> compute_freq('ACGCGGCTCTGAAA', 2)
        [2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]
    """
    freq_arr = [0] * (4**k)
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        num = pattern_to_number(kmer)
        freq_arr[num] += 1
    return freq_arr


# BA01K: compute frequency
def compute_freq_dict(Text, k):
    """
    (str,int) -> {int:int}
    >>> compute_freq_dict('ACGCGGCTCTGAAA', 2)
        {0:2,1:1,2:0,3:0,4:0,5:0,6:2,7:2,8:1,9:2,10:1,11:0,12:0,13:1,14:1,15:0}
    """
    FrequencyArray = {}    # dictionary, not list
    for i in range(0, 4**k):
        FrequencyArray[i] = 0
    for i in range(0, len(Text)-k+1):
        Pattern = Text[i:i+k]
        j = pattern_to_number(Pattern)
        FrequencyArray[j] += 1
    return FrequencyArray


# BA01B: frequent words
def faster_freq_words(text, k):
    """
    (str,int) -> {str}
    >>> faster_freq_words('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
        {'GCAT','CATG'}
    """
    freq_patterns = set()
    freq_arr = compute_freq(text, k)
    max_freq = max(freq_arr)
    for i in range(4**k):
        if freq_arr[i] == max_freq:
            kmer = number_to_pattern(i, k)
            freq_patterns.add(kmer)
    return freq_patterns


#BA01B: frequent words, faster
def faster_freq_words_dict(text, k):
    """
    (str,int) -> {str}
    >>> faster_freq_words_dict('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
        {'GCAT','CATG'}
    """
    freq_patterns = set()
    freq_dict = compute_freq_dict(text, k)
    max_freq = max(freq_dict.values())
    for i in range(4**k):
        if freq_dict[i] == max_freq:
            kmer = number_to_pattern(i, k)
            freq_patterns.add(kmer)
    return freq_patterns


# main
def main():
    # read lines
    f = open('/home/wsl/rosalind/data/ba01b.txt', 'r')
    lines = f.readlines()
    dna = lines[0]
    k = int(lines[1])
    f.close()

    # print answer
    start_time = time.time()
    print(faster_freq_words(dna, k))
    print("--- %s seconds ---" % (time.time() - start_time))  # --- 19.79618000984192 seconds ---

    start_time = time.time()
    #print(faster_freq_words_dict(dna, k))
    print("--- %s seconds ---" % (time.time() - start_time))  # --- 21.159732341766357 seconds ---


# execute
if __name__ == "__main__":
    main()
