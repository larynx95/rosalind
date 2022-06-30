/*
Rosalind: BA3B
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern[1], ... , Pattern[n]
such that the last k - 1 symbols of Pattern[i] are equal
to the first k - 1 symbols of Pattern[i+1] for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> PREV: Generate the k-mer Composition of a String (BA3A)
      -> HERE: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * NEXT: Construct OverapGraph (BA3C)

Plan 1.
  * sample dataset
    suffix --> prefix ...
    TAATGCCATGGGATGTT
    TAA->AAT->ATG->TGC->GCC->CCA->CAT->ATG->TGG->GGG->GGA->GAT->ATG->TGT->GTT
*/

// BA03B: genome path
function genome_path(texts) {
    /**
     * [str] -> str
     * returns a string from k-mer strings
     */
    let result = texts[0];
    let overlap = texts[0].slice(1);
    for (let i = 1; i < texts.length; i++) {
        result += texts[i].slice(texts[i].length - 1);
    }
    return result;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03b.txt').toString().split("\n");

    const startTime = performance.now()
    const result = genome_path(lines);
    console.log(result);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()