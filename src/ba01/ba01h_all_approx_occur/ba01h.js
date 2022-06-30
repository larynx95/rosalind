/*
Rosalind: BA1H
Find All Approximate Occurrences of a Pattern in a String

We say that a k-mer Pattern appears as a substring of Text with at most d mismatches
if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern,
i.e., HammingDistance(Pattern, Pattern') ≤ d.
Our observation that a DnaA box may appear with slight variations
leads to the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return:
All starting positions where Pattern appears
as a substring of Text with at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3

Sample Output
6 7 26 27 78

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
      -> PREV: hamming distance (BA1G)
      -> HERE: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> frequency array
          -> NEXT: number to pattern (BA1M), pattern to number (BA1L)

Plan 1.
- pseudocode:
    APPROXIMATEPATTERNCOUNT(Text, Pattern, d)
      count 0
      for i <- 0 to |Text| - |Pattern|
        Pattern' <- Text(i, |Pattern|)
        if HAMMINGDISTANCE(Pattern, Pattern') >= d
          count count + 1
      return count
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

// BA01H: approximate pattern count
function approx_pattern_count(text, pattern, d) {
    /**
     * (str,str,int) -> [int]
     * returns a list of starting indices
     */
    let indices = [];
    const k = pattern.length;
    for (let i = 0; i < text.length - pattern.length + 1; i++) {
        let pat = text.slice(i, i + k);
        if (hamming_distance(pat, pattern) <= d) {
            indices.push(i);
        }
    }
    return indices;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01h.txt').toString().split("\n");
    const pattern = lines[0];
    const text = lines[1];
    const d = parseInt(lines[2]);

    const startTime = performance.now();
    let indices = approx_pattern_count(text, pattern, d);
    answer = "";
    for (const idx of indices)
        answer += idx + " ";
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
