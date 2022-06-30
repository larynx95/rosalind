/*
Rosalind: BA1C
Find the Reverse Complement of a String

In DNA strings, symbols 'A' and 'T' are complements of each other,
as are 'C' and 'G'.
Given a nucleotide p, we denote its complementary nucleotide as p.
The reverse complement of a DNA string Pattern = p1…pn is the string rPattern = pn … p1
formed by taking the complement of each nucleotide in Pattern,
then reversing the resulting string.

For example, the reverse complement of Pattern = "GTCA" is rPattern = "TGAC".

Reverse Complement Problem
Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: rPattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT

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
      -> Find the Most Frequent Words with Mismatches in a String (BA1I)
      -> NEXT: Find Frequent Words with Mismatches and Reverse Complements (BA1)
        -> HERE: Find the Reverse Complement of a String (BA1C)

Plan 1.
- human DNA is double stranded helical structure, and has direction

  5-AAAACCCGGT-3
  3-  ...     -5     I should read this strand by 5 to 3 direction.

- steps:
  a. reverse given DNA strand
  b. get complement nucleotide for each character, that's all

═════════════════════════════════════════════════

References:
- How can I reverse an array in JavaScript without using libraries?
  https://stackoverflow.com/questions/10168034/how-can-i-reverse-an-array-in-javascript-without-using-libraries
- How do you reverse a string in-place in JavaScript?
  https://stackoverflow.com/questions/958908/how-do-you-reverse-a-string-in-place-in-javascript
*/

// BA01C: reverse complement
function reverse_complement(pattern) {
    /**
     * str -> str
     * returns a reverse complement (BA1C)
     * reverse_complement("ACGGTT") == AACCGT
     */
    dic = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' };
    pattern = [...pattern].reverse().join("");
    answer = "";
    for (const chr of pattern) {
        answer += dic[chr];
    }
    return answer;
}

// main
function main() {
    const fs = require('fs');
    let lines = fs.readFileSync('/home/wsl/rosalind/data/ba01c.txt').toString().split("\n");
    text = lines[0]

    let startTime = performance.now()
    console.log(reverse_complement(text));
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
