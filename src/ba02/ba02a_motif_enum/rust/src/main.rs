/*
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

  ╔═════════════════════════════════════════════════════════════════════════════════════╗
  ║  MOTIFENUMERATION(Dna, k, d)                                                        ║
  ║      Patterns <- an empty set                                                       ║
  ║      for each k-mer Pattern in Dna                                                  ║
  ║          for each k-mer Pattern' differing from Pattern by at most d mismatches     ║
  ║              if Pattern' appears in each string from Dna with at most d mismatches  ║
  ║                  add Pattern' to Patterns                                           ║
  ║      remove duplicates from Patterns                                                ║
  ║      return Patterns                                                                ║
  ╚═════════════════════════════════════════════════════════════════════════════════════╝

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
- What exactly is the meaning of this exercise?
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
  - kmers -> neighbors -> set intersection
  - get all neighberhood from all k-mers in the first dna string
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

References:
- Merge two HashMaps in Rust
  https://stackoverflow.com/questions/27244465/merge-two-hashmaps-in-rust
- How can I insert all values of one HashSet into another HashSet?
  https://stackoverflow.com/questions/56148586/how-can-i-insert-all-values-of-one-hashset-into-another-hashset
- intersection of two HashSets
  std::collections::hash_set::Intersection
  https://doc.rust-lang.org/std/collections/hash_set/struct.Intersection.html
  ; create a new HashSet
*/

// #!/usr/bin/env rust
use std::collections::HashSet;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01G: hamming distance
fn hamming_distance(astr: &str, bstr: &str) -> usize {
    astr.char_indices()
        .filter(|&(i, ch)| ch != bstr.chars().nth(i).unwrap())
        .map(|(_, ch)| ch)
        .count()
}

// BA01N: recursive neighbors
fn neighbors(pattern: &str, d: usize) -> HashSet<String> {
    if d == 0 {
        return HashSet::from([pattern.to_string()]);
    }
    if pattern.len() == 1 {
        return HashSet::from([
            "A".to_string(),
            "C".to_string(),
            "G".to_string(),
            "T".to_string(),
        ]);
    }
    let mut neighborhood: HashSet<String> = HashSet::new();
    let suffix_neighbors = neighbors(&pattern[1..], d);
    for text in suffix_neighbors {
        if hamming_distance(&pattern[1..], &text) < d {
            for nuc in "ACGT".chars() {
                neighborhood.insert(nuc.to_string() + &text);
            }
        } else {
            neighborhood.insert(pattern.chars().nth(0).unwrap().to_string() + &text);
        }
    }
    return neighborhood;
}

// BA02A: MotifEnumeration
fn motif_enumeration(dnas: &[String], k: usize, d: usize) -> HashSet<String> {
    let mut candidates: HashSet<String> = HashSet::new();
    // get a set of candidates from the first dna string
    for i in 0..dnas[0].len() - k + 1 {
        let kmer = &dnas[0][i..i + k]; // k-mer
        let neighbors = neighbors(kmer, d);
        candidates.extend(neighbors);
    }
    // get all neighbors from each dna string, then set intersection
    for dna in &dnas[1..] {
        // get all neighbors (neighborhood) from each dna string
        let mut neighborhood: HashSet<String> = HashSet::new();
        for i in 0..dna.len() - k + 1 {
            let kmer = &dna[i..i + k];
            let neighbors = neighbors(kmer, d);
            neighborhood.extend(neighbors);
        }
        // set inetersection, TODO: This part is not easy!
        // use 'retain' function
        candidates = candidates
            .iter()
            .filter(|n| neighborhood.contains(n.as_str()))
            .map(|str| str.to_string())
            .collect::<HashSet<String>>();
    }
    return candidates;
}

fn main() {
    // read text file, assign data to variables
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02a.txt");
    let split = lines[0]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let d = split[1];
    let dnas = &lines[1..];
    // dbg!(k, d, dnas);

    // check time
    let start = Instant::now();
    for motif in motif_enumeration(dnas, k, d) {
        print!("{} ", motif)
    }
    println!();
    println!("{:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
