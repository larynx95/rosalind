/*
Rosalind: BA1M
Implement n2p

Convert an integer to its corresponding DNA string.

Given: Integers index and k.
Return: n2p(index, k).

Sample Dataset
45
4

Sample Output
AGTC

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
  ║                                                                    ║
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
  ║                                                          ║
  ╚══════════════════════════════════════════════════════════╝
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

// BA01M: number to pattern, recursive
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

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01m.txt').toString().split("\n");
    const number = parseInt(lines[0]);
    const k = parseInt(lines[1]);

    let startTime = performance.now()
    console.log(number_to_pattern(number, k));
    console.log(`${performance.now() - startTime} milliseconds`)

    startTime = performance.now()
    console.log(number_to_pattern_rec(number, k));
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
