/*
Rosalind: BA1G
Compute the Hamming Distance Between Two Strings

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
For example, CGAAT and CGGAC have two mismatches. The number of mismatches
between strings p and q is called the Hamming distance between these strings and
is denoted HammingDistance(p, q).

Hamming Distance Problem
Compute the Hamming distance between two DNA strings.

Given: Two DNA strings.

Return: An integer value representing the Hamming distance.

Sample Dataset
GGGCCGTTGGT
GGACCGTTGAC

Sample Output
3

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
      -> PREV: Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> HERE: hamming distance (BA1G)
      -> NEXT: Find All Approximate Occurrences of a Pattern in a String (BA1H)
*/

// BA01G: hammind distance
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

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01g.txt').toString().split("\n");
    const astr = lines[0];
    const bstr = lines[1];

    const startTime = performance.now();
    console.log(hamming_distance(astr, bstr));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
