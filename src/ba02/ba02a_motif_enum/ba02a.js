/*
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

╔════════════════════════════════════════════════════════════════════════════════════╗
║ MOTIFENUMERATION(Dna, k, d)                                                        ║
║     Patterns <- an empty set                                                       ║
║     for each k-mer Pattern in Dna                                                  ║
║         for each k-mer Pattern' differing from Pattern by at most d mismatches     ║
║             if Pattern' appears in each string from Dna with at most d mismatches  ║
║                 add Pattern' to Patterns                                           ║
║     remove duplicates from Patterns                                                ║
║     return Patterns                                                                ║
╚════════════════════════════════════════════════════════════════════════════════════╝

Implanted Motif Problem
Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.

Given: Integers k and d, followed by a collection of strings Dna.

Return: All (k, d)-motifs in Dna.

Sample Dataset
3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output
ATA ATT GTT TTT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    * find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> HERE: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference?
      -> NEXT: "Median string problem" (BA2B)

Info:
                          T  C  G  G  G  G  g  T  T  T  t  t
                          c  C  G  G  t  G  A  c  T  T  a  C
                          a  C  G  G  G  G  A  T  T  T  t  C
                          T  t  G  G  G  G  A  c  T  T  t  t
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C
                          T  t  G  G  G  G  A  c  T  T  C  C
                          T  C  G  G  G  G  A  T  T  c  a  t
                          T  C  G  G  G  G  A  T  T  c  C  t
                          T  a  G  G  G  G  A  a  c  T  a  C
                          T  C  G  G  G  t  A  T  a  a  C  C
    SCORE(Motifs)         3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

Thoughts:
- What exactly is the meaning of this exercise? What are the problems?
  - finding motifs from multiple DNA string (e.g. 15-mer with mutations from each of ten DNA string )
  - Can this exercise be solved by functions in BA1B or BA1I? No.
    cuz ...
       i.   too slow
       ii.  mutations --> BA1B (x)
       iii. ten randomly generated DNA strings
            --> BA1I (x), solution of BA1I is only for single text
            --> concatenate ten DNA strings is indadequate (incorrect model of biological problem)

Terminology:
- motifs (k-mers)
  ; candiates for regulatory motif (or transcription factor binding site)
  ; a collection of k-mers, one from each string in DNA,
    minimizing SCORE(Motifs) among all possible choices of k-mers
- d(Pattern, text)
  ┌──────────────────────────────────────────────────────────────────────┐
  │  d(Pattern, Text) = min          HammingDistance(Pattern, Pattern')  │
  │                  all k-mers Pattern' in Text                         │
  └──────────────────────────────────────────────────────────────────────┘
- d(Pattern, Motifs)
  ┌───────────────────────────────────────────┐
  │                  t                        │
  │  d(Pattern, Dna) = Σ  d(Pattern, Dna_i)   │
  │                  i=0                      │
  └───────────────────────────────────────────┘

Plan 1.
  * brute-force
  * steps:
    a. candidates: all 4^k k-mer patterns as in frequency array (instead of k-mers from a string)
    b. find all (k,d)-mismatched k-mers from the first DNA string
    c. for the rest of DNA strings, do the same thing as in step (b), and intersect two sets

Plan 2.
  * using BA1N neighbors function
  * sample dataset
    dnas     k=3    d=1                                      motifs
    -----------------------------------------------------------------
    ATTTGGC  ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT  ATA ATT GTT TTT
             TTT    TTA TTC TTG TAT TCT TGT ATT CTT GTT TTT
             TTG    TTA TTC TAG TCG TGG ATG CTG GTG TTG TTT
             TGG    TGA TGC TAG TCG AGG CGG GGG TGG TTG TGT
             GGC    GGA GAC GCC AGC CGC GGC TGC GTC GGG GGT
    TGCCTTA  TGC    TGA TAC TCC AGC CGC GGC TGC TTC TGG TGT
             GCC    GCA GAC ACC CCC GCC TCC GGC GTC GCG GCT
             CCT    CCA CCC CCG CAT ACT CCT GCT TCT CGT CTT
             CTT    CTA CTC CTG CAT CCT CGT ATT CTT GTT TTT
             TTA    TAA TCA TGA ATA CTA GTA TTA TTC TTG TTT
    CGGTATC  CGG    CGA CGC CAG CCG AGG CGG GGG TGG CTG CGT
             GGT    GGA GGC GGG GAT GCT AGT CGT GGT TGT GTT
             GTA    GAA GCA GGA ATA CTA GTA TTA GTC GTG GTT
             TAT    TAA TAC TAG AAT CAT GAT TAT TCT TGT TTT
             ATC    ATA AAC ACC AGC ATC CTC GTC TTC ATG ATT
    GAAAATT  GAA    AAA CAA GAA TAA GCA GGA GTA GAC GAG GAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAT    AAA AAC AAG AAT CAT GAT TAT ACT AGT ATT
             ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT

═════════════════════════════════════════════════

References:
- Javascript set operations
- higher order function
*/

// BA01G: hamming distance
function hamming_distance(astr, bstr) {
    /**
     * (str,str) -> int
     * returns a Hamming distance of two strings (BA1G)
     */
    let distance = 0;
    for (let i = 0; i < astr.length; i++)
        if (astr[i] != bstr[i]) distance++;
    return distance;
}

// BA01N: neighbors
function neighbors(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors (BA1N)
     */
    // Part A. base case
    if (d == 0) return new Set([pattern]);
    if (pattern.length == 0) return new Set();
    if (pattern.length == 1) return new Set(['A', 'C', 'G', 'T']);
    // Part B. recursive case
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

// BA02A: motif enumeration
function motif_enum(dnas, k, d) {
    /**
     * ([str],int,int) -> {str}
     * returns a set of all (k, d)-motifs (BA2A)
     * TODO: Is this function efficient?
     */
    // get all d-neighbors from each dna string
    let ls_candidates = [];
    for (const dna of dnas) {
        let candidates = new Set([]);
        for (let i = 0; i < dna.length - k + 1; i++) {
            const nbs = neighbors(dna.slice(i, i + k), d);
            candidates = new Set([...candidates, ...nbs]);
        }
        ls_candidates.push(candidates);
    }
    // set intersection
    const answer = new Set(ls_candidates.reduce((a, b) => [...a].filter(x => b.has(x))));
    return answer;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02a.txt').toString().split("\n");
    const [k, d] = lines[0].split(/\s/).map(Number);
    const dnas = lines.slice(1);

    const startTime = performance.now();
    const motifs = motif_enum(dnas, k, d);
    answer = "";
    for (const motif of motifs) {
        answer += motif + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
