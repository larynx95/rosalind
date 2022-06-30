/*
Rosalind: BA1E
Find Patterns Forming Clumps in a String

Given integers L and t,
a string Pattern forms an (L, t)-clump inside a (larger) string Genome
if there is an interval of Genome of length L
in which Pattern appears at least t times.
For example, TGCA forms a (25,3)-clump in the following Genome:
gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

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
  ║      if FREQUENCYARRAY(index) >= t                        ║
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

    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT AAGCAAAGGTGGGC   fst  lst
    3  0  2  0  1  0  0  0  0  1  1  1  0  0  1  0  AAGCAAAGGTG
   -2  0  2  0  1  0  0  0  0  1  2+ 1  0  0  1  0   AGCAAAGGTGG     AA   GG
    2  0 -1  0  1  0  0  0  0  1  3+ 1  0  0  1  0    GCAAAGGTGGG    AG   GG
    2  0  1  0  1  0  0  0  0 -1+ 3  1  0  0  1  0     CAAAGGTGGGC   GC   GC
    check the frequency before subtraction!

═════════════════════════════════════════════════

References:
- Does JavaScript support array/list comprehensions like Python?
  https://stackoverflow.com/questions/31353213/does-javascript-support-array-list-comprehensions-like-python
- How to split a string by white space or comma?
  https://stackoverflow.com/questions/10346722/how-to-split-a-string-by-white-space-or-comma?rq=1
- How to convert a string of numbers to an array of numbers?
  https://stackoverflow.com/questions/15677869/how-to-convert-a-string-of-numbers-to-an-array-of-numbers
  - var b = a.split(',').map(Number);
  - parseInt(item, 10);
  - ["1", "2", "3", "4"].map(i=>Number(i))
*/

// BA01M: number ot pattern
function number_to_pattern(num, k) {
    /**
     * (int,int) -> str
     * returns a pattern from a number
     * number_to_pattern(15, 2) == "TT"
     */
    pattern = ""
    for (let i = 0; i < k; i++) {
        temp = Math.floor(num / Math.pow(4, (k - i - 1)));
        if (temp == 0) pattern += 'A';
        else if (temp == 1) pattern += 'C';
        else if (temp == 2) pattern += 'G';
        else if (temp == 3) pattern += 'T';
        num -= temp * Math.pow(4, k - i - 1);
    }
    return pattern
}

// BA01L: pattern to number
function pattern_to_number(pattern) {
    /**
     * str -> int
     * returns a number from a pattern
     * pattern_to_number("TT") == 15
     */
    const dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
    let number = 0;
    for (let i = 0; i < pattern.length; i++) {
        let pw = pattern.length - i - 1;
        let nuc = 0;
        let chr = pattern[i];
        if (chr == 'A') nuc = 0;
        else if (chr == 'C') nuc = 1;
        else if (chr == 'G') nuc = 2;
        else if (chr == 'T') nuc = 3;
        let pw_val = Math.pow(4, pw);
        number += pw_val * nuc;
    }
    return number;
}

// BA01E
function find_freq_patterns(text, k, t) {
    /**
     * (str,int,int) -> {str}
     * returns a frequency array (modified "compute_freq_dic" function)
     */
    let freq_dic = {};
    let freq_kmers = new Set();
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        if (!(num in freq_dic)) freq_dic[num] = 1;
        else {
            freq_dic[num] += 1;
            if (freq_dic[num] >= t) freq_kmers.add(number_to_pattern(num, k));
        }
    }
    return freq_kmers
}

// BA01E: clume finding
function clump_finding(genome, k, l, t) {
    /**
     * (str,int,int,int) -> {str}
     * returns k-mers forming (L, t)-clumps
     */
    let freq_kmers = new Set();
    for (let i = 0; i < genome.length - l + 1; i++) {
        let text = genome.slice(i, i + l);
        find_freq_patterns(text, k, t).forEach(freq_kmers.add, freq_kmers);
    }
    return freq_kmers;
}

// BA01K: compute frequency
function compute_freq_dic(text, k) {
    /**
     * (str,int) -> {int:int}
     * returns a frequency array
     */
    let freq_dic = {};
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        if (!(num in freq_dic)) freq_dic[num] = 1;
        else ++freq_dic[num];
    }
    return freq_dic;
}

// BA01E: clump finding better
function clump_finding_better(genome, k, l, t) {
    /**
     * (str,int,int,int) -> {str}
     * returns k-mers forming (L, t)-clumps
     */
    // construct frequency array
    let freq_kmers = new Set();
    let freq_dic = compute_freq_dic(genome.slice(0, l), k);
    for (let i = 1; i < genome.length - l + 1; i++) {
        // add the last k-mer in the new clump, before check frequency (The order matters.)
        let lst = pattern_to_number(genome.slice(i + l - k, i + l));
        let fst = pattern_to_number(genome.slice(i - 1, i + k - 1));
        if (!(lst in freq_dic)) freq_dic[lst] = 1
        else ++freq_dic[lst];
        if (freq_dic[fst] >= t) {
            freq_kmers.add(number_to_pattern(fst, k));
        }
        // remove the first k-mer in the previous clump, after check frequency (The order matters.)
        if (freq_dic[lst] >= t) {
            freq_kmers.add(number_to_pattern(lst, k));
        }
        --freq_dic[fst];
    }
    return freq_kmers;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01e.txt').toString().split("\n");
    const genome = lines[0];
    const [k, l, t] = lines[1].split(/\s+/).map(Number);  // useful expression

    // clump_finding
    let startTime = performance.now();
    let patterns = clump_finding(genome, k, l, t);
    let answer = "";
    for (const pat of patterns) {
        answer += pat + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);  // 5487.461899999529 milliseconds

    // clump_finding_better (much faster!)
    startTime = performance.now();
    patterns = clump_finding_better(genome, k, l, t);
    answer = "";
    for (const pat of patterns) {
        answer += pat + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);  // 33.400600001215935 milliseconds
}

// execute main function
main()
