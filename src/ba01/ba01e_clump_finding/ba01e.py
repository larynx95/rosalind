"""
Rosalind: BA1E
Find Patterns Forming Clumps in a String

Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger)
string Genome if there is an interval of Genome of length L in which Pattern
appears at least t times. For example, TGCA forms a (25,3)-clump in the
following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem
Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset
CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4

Sample Output
CGACA GAAGA AATGT

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
      -> PREV: find all occurrence of a pattern in a string (BA1D)
      -> HERE: Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> NEXT: hamming distance (BA1G)

Plan 1.
- brute-force

Plan 2.
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  CLUMPFINDING(Genome, k, t, L)                                        ║
  ║    FrequentPatterns <- an empty set                                   ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLUMP(i) <- 0                                                    ║
  ║    for i <- 0 to |Genome| - L                                         ║
  ║      Text <- the string of length L starting at position i in Genome  ║
  ║      FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)                  ║
  ║      for index <- 0 to 4k - 1                                         ║
  ║        if FREQUENCYARRAY(index) >= t                                  ║
  ║          CLUMP(index) <- 1                                            ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLUMP(i) = 1                                                  ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                               ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

Plan 3.
  ╔═══════════════════════════════════════════════════════════╗
  ║  CLUMPFINDINGBetter(Genome, k, t, L)                      ║
  ║    FrequentPatterns <- an empty set                       ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      CLUMP(i) <- 0                                        ║
  ║    Text <- Genome(0, L)                                   ║
  ║    FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)        ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if FREQUENCYARRAY(i) - t                             ║
  ║        CLUMP(i) <- 1                                      ║
  ║    for i <- 1 to |Genome| - L                             ║
  ║      FirstPattern <- Genome(i - 1, k)                     ║
  ║      index <- PATTERNTONUMBER(FirstPattern)               ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) - 1   ║
  ║      LastPattern <- Genome(i + L - k, k)                  ║
  ║      index <- PATTERNTONUMBER(LastPattern)                ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) + 1   ║
  ║      if FREQUENCYARRAY(index) - t                         ║
  ║        CLUMP(index) <- 1                                  ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if CLUMP(i) = 1                                      ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                   ║
  ║        add Pattern to the set FrequentPatterns            ║
  ║    return FrequentPatterns                                ║
  ╚═══════════════════════════════════════════════════════════╝

1. find the generalized rule of frequency array

  "AAGCAAAGGTGGGC"  len=14
   vv-----------                                    first ... last
  "AAGCAAAGGTGGG"   len=13, starting from index 0   AA    ...  -
   "AGCAAAGGTGGGC"  len=13, starting from index 1   -     ...  GC
    -----------^^
    common part!

  (1) compute_freq("AAGCAAAGGTGGG", 2)
                    ^^  <--- minus 1
  (2) compute_freq( "AGCAAAGGTGGGC", 2)
                                ^^  <--- plus 1

      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
      3  0  2  0  1  0  0  0  0  1  3  1  0  0  1  0   --- (1)
     -2  0  2  0  1  0  0  0  0  2+ 3  1  0  0  1  0   --- (2)

2. If we know a frequency array of the first clump,
   we can get the frequency array of the whole genome.

    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT AAGCAAAGGTGGGC
    3  0  2  0  1  0  0  0  0  1  1  1  0  0  1  0  AAGCAAAGGTG
   -2  0  2  0  1  0  0  0  0  1  2+ 1  0  0  1  0   AGCAAAGGTGG
    2  0 -1  0  1  0  0  0  0  1  3+ 1  0  0  1  0    GCAAAGGTGGG
    2  0  1  0  1  0  0  0  0 -1+ 3  1  0  0  1  0     CAAAGGTGGGC
"""

import time


#################################################
# brute-force algorithm
#################################################

# BA01D: find all occurrences
def find_all_occur(pattern, dna):
    """
    (str,str) -> [int]
    >>> all_occurrence('ATAT','GATATATGCATATACTT')
        [1,3,9]
    """
    indices = []
    for i in range(0, len(dna)-len(pattern)+1):
        if pattern == dna[i:i+len(pattern)]:
            indices.append(i)
    return indices


# BA01E: clump finding
def clump_finding_brute(genome, k, L, t):
    """
    (str,int,int,int) -> {str}
    >>> genome = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
    >>> clump_finding_brute(genome,5,75,4)
        {'GAAGA','CGACA','AATGT'}
    """
    patterns = set()
    # outer loop: for all possible clumps
    for i in range(0, len(genome) - L + 1):
        clump = genome[i:i+L]
        # inner loop: for all kmers
        for j in range(0, len(clump) - k + 1):
            kmer = clump[j:j+k]
            freq = len(find_all_occur(kmer, clump))
            if freq >= t:
                patterns.add(kmer)
    return patterns


#################################################
# clump finding by frequency array
#################################################

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
    (str,int) -> {int:int}
    >>> compute_freq('ACGCGGCTCTGAAA', 2)
        {0:2,1:1,2:0,3:0,4:0,5:0,6:2,7:2,8:1,9:2,10:1,11:0,12:0,13:1,14:1,15:0}
    """
    freq_arr = {}    # dictionary, not list
    for i in range(0, 4**k):
        freq_arr[i] = 0
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        freq_arr[j] += 1
    return freq_arr


# BA01E: clump finding
def clump_finding(genome,k,l,t):
    """
    (str,int,int,int) -> {str}
    >>> genome = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
    >>> clump_finding(genome,5,75,4)
        {'AATGT','GAAGA','CGACA'}
    """
    freq_patterns = set()
    clump = [0] * 4**k
    for i in range(len(genome)-l+1):
        text = genome[i:i+l]
        freq_arr = compute_freq(text,k)
        for i in range(4**k):
            if freq_arr[i] >= t:
                clump[i] = 1
    for i in range(4**k):
        if clump[i] == 1:
            pattern = number_to_pattern(i,k)
            freq_patterns.add(pattern)
    return freq_patterns


#################################################
# clump finding by better algorithm
#################################################

# BA01E: clump finding, better
def clump_finding_better(Genome, k, L, t):
    FrequentPatterns = []
    CLump = {}
    # initialize an empty "dictionary(fast)" with zero
    for i in range(0, 4**k):
        CLump[i] = 0

    # the first clump, and its frequency array
    Text = Genome[0:L+1]
    FrequencyArray = compute_freq(Text, k)

    # for the first clump,
    #   a. with frequency array for the first clump
    #   b. get patterns (as numbers) satisfying frequency condition
    for index in range(4**k):
        if FrequencyArray[index] >= t:
            CLump[index] = 1

    # TODO: mysterious part:
    # Key concept:
    #   If we know the frequency array of the first clump,
    #   we can get frequency arrays for the other clumps.
    #   The only difference between two consecutive strings is
    #       - the first part of the previous string
    #       - and the last part of the next string.
    # Advantages:
    #   In this approach,
    #   we just add/subtract frequencies of the the first/last kmers.
    #   That is, we can ignore (or skip) frequencies of kmers from middle part.
    for i in range(1, len(Genome)-L+1):
        FirstPattern = Genome[i-1:(i-1)+k]
        index = pattern_to_number(FirstPattern)
        FrequencyArray[index] -= 1
        LastPattern = (Genome[i + L - k:(i + L - k)+k])
        index = pattern_to_number(LastPattern)
        FrequencyArray[index] += 1
        if FrequencyArray[index] >= t:
            CLump[index] = 1

    for i in range(4**k):
        if CLump[i]==1:
            Pattern = number_to_pattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01e.txt', 'r')
    lines = f.readlines()
    genome = lines[0].strip()
    numbers = [int(line.strip()) for line in lines[1].split()]
    k = numbers[0]
    L = numbers[1]
    t = numbers[2]
    f.close()

    start_time = time.time()
    # print(clump_finding_brute(genome, k, L, t))                      # --- 458.49220633506775 seconds ---
    print(clump_finding_better(genome, k, L, t))
    print (' '.join(map(str, clump_finding_better(genome, k, L, t))))  # --- 1.8848645687103271 seconds ---
    print("--- %s seconds ---" % (time.time() - start_time))


# execute
if __name__ == "__main__":
    main()