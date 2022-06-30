/*
Rosalind: BA3A (difficulty: 1/5)
Generate the k-mer Composition of a String

Given a string Text, its k-mer composition Compositionk(Text)
is the collection of all k-mer substrings of Text
(including repeated k-mers).
For example,

Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}

Note that we have listed k-mers in lexicographic order
(i.e., how they would appear in a dictionary)
rather than in the order of their appearance in TATGGGGTGC.
We have done this because the correct ordering of the reads
is unknown when they are generated.

String Composition Problem
Generate the k-mer composition of a string.

Given: An integer k and a string Text.

Return: Compositionk(Text)
        (the k-mers can be provided in any order).

Sample Dataset
5
CAATCCAAC

Sample Output
AATCC
ATCCA
CAATC
CCAAC
TCCAA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> HERE: Generate the k-mer Composition of a String (BA3A)
      -> NEXT: String Reconstruction Problem with k-mer composition(BA3B)
*/

// BA03A: all k-mers
function all_kmers_clang_style(text, k) {
    /**
     * (str,int) -> [str]
     * returns a list of all k-mers (c/cpp style)
     */
    let kmers = [];
    for (let i = 0; i < text.length; i++) {
        kmers.push(text.slice(i, i + k));
    }
    return kmers;
}

// BA03A: all k-mers
function* all_kmers_gen(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length; i++) yield text.slice(i, i + k);
}

// BA03: all k-mers
function all_kmers(text, k) {
    /**
     * (str,int) -> [str]
     * returns a list of all k-mers
     */
    return [...Array(text.length - k + 1)].map((_, i) => e = text.slice(i, i + k));
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03a.txt').toString().split("\n");
    const k = parseInt(lines[0]);
    const text = lines[1];

    const startTime = performance.now();
    const kmers = all_kmers(text, k);
    for (const kmer of kmers) {
        console.log(kmer);
    }
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()