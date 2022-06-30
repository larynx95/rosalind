/*
Rosalind: BA1D
Find All Occurrences of a Pattern in a String

In this problem, we ask a simple question:
how many times can one string occur as a substring of another?
Recall from “Find the Most Frequent Words in a String”
that different occurrences of a substring can overlap with each other.
For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
Find all occurrences of a pattern in a string.

Given:
Strings Pattern and Genome.

Return:
All starting positions in Genome where Pattern appears as a substring.
Use 0-based indexing.

Sample Dataset
ATAT
GATATATGCATATACTT

Sample Output
1 3 9

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
      -> HERE: find all occurrence of a pattern in a string (BA1D)
      -> NEXT: Clump finidng problem (BA1E)

═════════════════════════════════════════════════

Refrences:
- What is the correct way to check for string equality in JavaScript?
  https://stackoverflow.com/questions/3586775/what-is-the-correct-way-to-check-for-string-equality-in-javascript
*/

// BA01D: all occurrences
function all_occurrences(pattern, genome) {
    /**
     * (str,str) -> [int]
     * returns a list of indices of pattern in a genome
     */
    let indices = [];
    let k = pattern.length;
    for (let i = 0; i < genome.length; i++) {
        if (genome.slice(i, i + k) === pattern) indices.push(i);
    }
    return indices
}

// main
function main() {
    const fs = require('fs');
    let lines = fs.readFileSync('/home/wsl/rosalind/data/ba01d.txt').toString().split("\n");
    pattern = lines[0]
    genome = lines[1]

    let startTime = performance.now()
    let ls = all_occurrences(pattern, genome);
    answer = "";
    for (const idx of ls) {
        answer += idx + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
