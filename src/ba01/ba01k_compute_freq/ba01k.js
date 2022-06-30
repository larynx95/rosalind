/*
Rosalind: BA1K
Generate the Frequency Array of a String

Given an integer k, we define the frequency array of a string Text as an array
of length 4k, where the i-th element of the array holds the number of times that
the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

kmer      AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
index      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
frequency  3   0   2   0   1   0   0   0   0   1   3   1   0   0   1   0

Computing a Frequency Array
Generate the frequency array of a DNA string.

Given: A DNA string Text and an integer k.

Return: The frequency array of k-mers in Text.

Sample Dataset
ACGCGGCTCTGAAA
2

Sample Output
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

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
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - by frequency array
        -> HERE: frequency array
          -> PREV: number to pattern (BA1M), pattern to number (BA1L)
        -> NEXT: neigbors (BA1N)

Plan 1.
- pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  COMPUTINGFREQUENCIES(Text, k)                    ║
  ║    for i <- 0 to 4^k- 1                           ║
  ║      FREQUENCYARRAY(i) <- 0                       ║
  ║    for i <- 0 to |Text| - k                       ║
  ║      Pattern <- Text(i, k)                        ║
  ║      j <- PATTERNTONUMBER(Pattern)                ║
  ║      FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1   ║
  ║    return FREQUENCYARRAY                          ║
  ╚═══════════════════════════════════════════════════╝

*/

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

// BA01L: number to pattern, recursive
function number_to_pattern_rec(num, k) {
    /**
     * (int,int) -> str
     * returns a pattern from a number
     * number_to_pattern_rec(15, 2) == "TT"
     */
    const dic = { 0: "A", 1: "C", 2: "G", 3: "T" };
    if (k == 0) return ""
    else {
        temp = Math.floor(num / Math.pow(4, (k - 1)));
        num -= temp * Math.pow(4, k - 1);
        return dic[temp] + number_to_pattern(num, k - 1)
    }
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

// BA01L: pattern to number, recursive
function pattern_to_number_rec(pattern) {
    /**
     * str -> int
     * returns a number from a pattern
     * pattern_to_number_rec("TT") == 15
     */
    const dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
    if (pattern.length == 0) {
        return 0;
    } else {
        return dic[pattern[0]] * Math.pow(4, pattern.length - 1) + pattern_to_number(pattern.slice(1));
    }
}

// BA01K: compute frequency
function compute_freq(text, k) {
    /**
     * (str,int) -> [int]
     * returns a frequency array
     * compute_freq('ACGCGGCTCTGAAA', 2)
     * == [2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]
     */
    let freq_arr = new Array(Math.pow(4, k));  // TODO: find the fastest way
    for (let i = 0; i < Math.pow(4, k); i++) freq_arr[i] = 0;
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        freq_arr[num] += 1;
    }
    return freq_arr;
}

// BA01K: compute frequency
function compute_freq_dic(text, k) {
    /**
     * (str,int) -> ({int:int},int)
     * returns a frequency array
     */
    let freq_dic = {};
    let max_cnt = 0;
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        if (!(num in freq_dic)) freq_dic[num] = 1;
        else {
            freq_dic[num] += 1;
            if (freq_dic[num] > max_cnt) max_cnt = freq_dic[num];
        }
    }
    return [freq_dic, max_cnt];
}

// main
function main() {
    const fs = require('fs');
    let lines = fs.readFileSync('/home/wsl/rosalind/data/ba01k.txt').toString().split("\n");
    const pattern = lines[0]
    const k = parseInt(lines[1]);

    let startTime = performance.now()
    let answer = "";
    for (const elem of compute_freq(pattern, k)) {
        answer += elem + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
