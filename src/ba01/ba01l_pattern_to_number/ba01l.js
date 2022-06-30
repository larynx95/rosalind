/*
Rosalind: BA1L
Implement PatternToNumber

Implement PatternToNumber
Convert a DNA string to a number.

Given: A DNA string Pattern.
Return: PatternToNumber(Pattern).

Sample Dataset
AGT

Sample Output
11

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
      -> PREV: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> NEXT: frequency array
          -> HERE: number to pattern (BA1M), pattern to number (BA1L)

Plan 1.
- pseudocode:
  ╔════════════════════════════════════════════════════════════════════╗
  ║  PATTERNTONUMBER(Pattern)                                          ║
  ║      if Pattern contains no symbols                                ║
  ║          return 0                                                  ║
  ║      symbol LASTSYMBOL(Pattern)                                    ║
  ║      Prefix PREFIX(Pattern)                                        ║
  ║      return 4 · PATTERNTONUMBER(Prefix) + SYMBOLTONUMBER(symbol)   ║
  ╚════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NUMBERTOPATTERN(index , k)                              ║
  ║      if k = 1                                            ║
  ║          return NUMBERTOSYMBOL(index)                    ║
  ║      prefixIndex QUOTIENT(index, 4)                      ║
  ║      r REMAINDER(index, 4)                               ║
  ║      symbol NUMBERTOSYMBOL(r)                            ║
  ║      PrefixPattern NUMBERTOPATTERN(prefixIndex, k  1)    ║
  ║      return concatenation of PrefixPattern with symbol   ║
  ╚══════════════════════════════════════════════════════════╝

*/

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

// BA01M: pattern to number
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

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01l.txt').toString().split("\n");
    const pattern = lines[0];

    let startTime = performance.now()
    console.log(pattern_to_number(pattern));
    console.log(`${performance.now() - startTime} milliseconds`)

    startTime = performance.now()
    console.log(pattern_to_number_rec(pattern));
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
