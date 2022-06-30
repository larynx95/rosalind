/*
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
  ║        index <- p2n(Pattern)                                          ║
  ║        CLOSE(index) <- 1                                              ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLOSE(i) = 1                                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        FREQUENCYARRAY(i) <- APPROXIMATEPATTERNCOUNT(Text, Pattern, d) ║
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
*/

// BA01G: hamming distance
function hamming_distance(astr, bstr) {
    /**
     * (str,str) -> int
     * returns a Hamming distance of two strings
     */
    let distance = 0;
    for (let i = 0; i < astr.length; i++)
        if (astr[i] != bstr[i]) distance++;
    return distance;
}

// BA01M: number to pattern
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

// BA01N: neighbors
function neighbors_iter(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    let neighborhood = new Set([pattern]);
    for (let i = 0; i < d; i++) {
        for (const pat of neighborhood) {
            neighborhood = new Set([...immediate_neighbors(pat), ...neighborhood]);
        }
    }
    return neighborhood;
}

// BA01N: neighbors
function neighbors(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    //
    if (d == 0) return new Set([pattern]);
    if (pattern.length == 0) return new Set();
    if (pattern.length == 1) return new Set(['A', 'C', 'G', 'T']);
    let neighborhood = new Set();
    let suffixneighbors = neighbors(pattern.slice(1), d);
    for (const text of suffixneighbors) {
        if (hamming_distance(text, pattern.slice(1)) < d) {
            for (const nuc of ['A', 'C', 'G', 'T']) {
                neighborhood.add(nuc + text);
            }
        } else {
            neighborhood.add(pattern[0] + text);
        }
    }
    return neighborhood;
}

// BA01I: frequent words mismatches
function freq_words_mismatches(text, k, d) {
    /**
     * (str,int,int) -> {str}
     * returns a set of frequent words with mismatches
     */
    let freq_patterns = new Set();
    let count = {};
    let max_cnt = -1;
    for (let i = 0; i < text.length - k + 1; i++) {
        let neighborhood = neighbors(text.slice(i, i + k), d);
        for (const word of neighborhood) {
            let idx = pattern_to_number(word);
            if (!(idx in count)) count[idx] = 1;
            else count[idx]++;
            if (count[idx] > max_cnt) max_cnt = count[idx];
        }
    }
    for (const num in count) {
        if (count[num] == max_cnt) {
            freq_patterns.add(number_to_pattern(num, k));
        }
    }
    return freq_patterns;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01i.txt').toString().split("\n");
    const text = lines[0]
    const [k, d] = [...lines[1].split(/\s/)].map(Number);

    const startTime = performance.now()
    const fwords = freq_words_mismatches(text, k, d);
    let answer = "";
    for (const word of fwords) {
        answer += word + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`)
    /*
    c++         0.4 sec
    go          0.102749 sec
    JavaScript  0.007 sec     <--- fastest! unexpected result!
    Python      0.56 sec

    Why does JavaScript appear to be 4 times faster than C++?
    https://stackoverflow.com/questions/17036059/why-does-javascript-appear-to-be-4-times-faster-than-c
    */
}

// execute main function
main()
