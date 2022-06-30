/*
Rosalind: BA1F
Find a Position in a Genome Minimizing the Skew

Define the skew of a DNA string Genome, denoted Skew(Genome),
as the difference between the total number of occurrences of 'G' and 'C' in Genome.
Let Prefixi(Genome) denote the prefix (i.e., initial substring) of Genome of length i.
For example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i
(from 0 to |Genome|).

Sample Dataset
CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

Sample Output
53 97

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> HERE: minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * NEXT: Most frequent words problem

Plan 1.
  * sample dataset
      C A T G G G C A T C G G C C A T A C G C C
     -1    +1+1+1-1    -1+1+1-1-1      -1+1-1-1
    0-1-1-1 0 1 2 1 1 1 0 1 2 1 0 0 0 0-1 0-1-2
*/

// BA01F: min skew
function min_skew(genome) {
    /**
     * str -> [int]
     * returns a list of indices in which skewness is minimum
     */
    let skewness = [0];
    const dic = { 'A': 0, 'C': -1, 'G': +1, 'T': 0 };
    let min_val = Infinity;
    for (const nuc of genome) {
        let cur_val = skewness[skewness.length - 1] + dic[nuc];
        skewness.push(cur_val);
        if (cur_val < min_val) min_val = cur_val;
    }
    let indices = [];
    for (let i = 0; i < skewness.length; i++) {
        if (min_val == skewness[i]) indices.push(i);
    }
    return indices;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01f.txt').toString().split("\n");
    const genome = lines[0]

    const startTime = performance.now()
    const indices = min_skew(genome);
    let answer = "";
    for (const idx of indices) {
        answer += idx + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
