/*
Rosalind: BA1N
Generate the d-Neighborhood of a String

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming
distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String
Find all the neighbors of a pattern.

Given: A DNA string Pattern and an integer d.
Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
ACG
1

Sample Output
CCG
TCG
GCG
AAG
ATG
AGG
ACA
ACC
ACT
ACG

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
      -> NEXT: Find the Most Frequent Words with Mismatches in a String (BA1I)
        -> PREV: frequency array
          -> number to pattern (BA1M), pattern to number (BA1L)
        -> HERE: neigbors (BA1N)

Plan 1.
  * immediate neighbors
    ACG, d=1
    A --> CCG, GCG, TCG
    C --> AAG, AGG, ATG
    G --> ACA, ACC, ACT
    ╔═══════════════════════════════════════════════════════════════════════╗
    ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
    ║    Neighborhood <- the set consisting of the single string Pattern    ║
    ║    for i = 1 to |Pattern|                                             ║
    ║      symbol <- i-th nucleotide of Pattern                             ║
    ║      for each nucleotide x different from symbol                      ║
    ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
    ║        add Neighbor to Neighborhood                                   ║
    ║    return Neighborhood                                                ║
    ╚═══════════════════════════════════════════════════════════════════════╝

  * iterative neighbors
    ╔══════════════════════════════════════════════════════════════╗
    ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
    ║    Neighborhood <- set consisting of single string Pattern   ║
    ║    for j = 1 to d                                            ║
    ║      for each string Pattern' in Neighborhood                ║
    ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
    ║        remove duplicates from Neighborhood                   ║
    ║    return Neighborhood                                       ║
    ╚══════════════════════════════════════════════════════════════╝

  * neighbors - recursive function
    (1) base case
      - if d == 0          => {pattern}
      - if |pattern| == 0  => {}
      - if |pattern| == 1  => {'A', 'C', 'G', 'T'}
    (2) recursive case
      - if neighbors(pattern[1:]) < d   => can prepend any nucleotide
      - if neighbors(pattern[1:]) >= d  => can change the first nucleotide
    ╔══════════════════════════════════════════════════════════╗
    ║  NEIGHBORS(Pattern, d)                                   ║
    ║    if d = 0                                              ║
    ║      return {Pattern}                                    ║
    ║    if |Pattern| = 1                                      ║
    ║      return {A, C, G, T}                                 ║
    ║    Neighborhood <- an empty set                          ║
    ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
    ║    for each string Text from SuffixNeighbors             ║
    ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
    ║        for each nucleotide x                             ║
    ║          add x + Text to Neighborhood                    ║
    ║      else                                                ║
    ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
    ║    return Neighborhood                                   ║
    ╚══════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- How to select all other values in an array except the ith element?
  https://stackoverflow.com/questions/15361189/how-to-select-all-other-values-in-an-array-except-the-ith-element
- do <something> N times (declarative syntax)
  https://stackoverflow.com/questions/10993824/do-something-n-times-declarative-syntax
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

// BA01N: immediate neighbors
function immediate_neighbors(pattern) {
    /**
     * str -> {str}
     * returns a set of immediated neighbors
     * immediate_neighbors("ACG")
     * == Set(10) {'ACG','CCG','GCG','TCG','AAG','AGG','ATG','ACA','ACC','ACT'}
     */
    neighborhood = new Set([pattern]);
    for (let i = 0; i < pattern.length; i++) {
        let symbol = pattern[i];
        for (const nuc of "ACGT".split('').filter(elem => elem != symbol)) {
            let neighbor = pattern.slice(0, i) + nuc + pattern.slice(i + 1);
            neighborhood.add(neighbor);
        }
    }
    return neighborhood;
}

// BA01N: iterative neighbors
function neighbors_iter(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    let neighborhood = new Set([pattern]);
    for (let i = 0; i < d; i++) {
        for (const pat of neighborhood) {
            neighborhood = new Set([...immediate_neighbors(pat), ...neighborhood]);
        }
    }
    return neighborhood;
}

// BA01N: recursive neighbors
function neighbors(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    //
    if (d == 0) return new Set([pattern]);
    if (pattern.length == 0) return new Set();
    if (pattern.length == 1) return new Set(['A', 'C', 'G', 'T']);
    let neighborhood = new Set();
    let suffixneighbors = neighbors(pattern.slice(1), d);
    for (const text of suffixneighbors) {
        if (hamming_distance(text, pattern.slice(1)) < d) {
            for (const nuc of ['A', 'C', 'G', 'T']) {
                neighborhood.add(nuc + text);
            }
        } else {
            neighborhood.add(pattern[0] + text);
        }
    }
    return neighborhood;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01n.txt').toString().split("\n");
    const pattern = lines[0];
    const d = parseInt(lines[1]);

    // print to screen
    const startTime = performance.now();
    const all_neighbors = neighbors(pattern, d);
    for (const neighbor of all_neighbors) {
        console.log(neighbor);
    }
    console.log(`${performance.now() - startTime} milliseconds`);

    // print to file
    // const file = fs.createWriteStream('../data/ba1n_output.txt');
    // all_neighbors.forEach((neighbor) => file.write(neighbor + '\n'));
    // file.end();
}

// execute main function
main()
