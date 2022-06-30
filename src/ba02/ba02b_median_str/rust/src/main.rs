/*
Rosalind: BA2B
Find a Median String

Given a k-mer Pattern and a longer string Text,
we use d(Pattern, Text) to denote the minimum Hamming distance
between Pattern and any k-mer in Text,

d(Pattern, Text) = min                 HammingDistance(Pattern, Pattern')
                  all k-mers Pattern' in Text

    example:
    d(GATTCTCA, gcaaaGACGCTGAccaa) = 3

Given a k-mer Pattern and a set of strings
Dna = {Dna_1, ... , Dna_t},
we define d(Pattern, Dna) as the sum of distances
between Pattern and all strings in Dna,

                    t
    d(Pattern, Dna) = Σ  d(Pattern, Dna_i)
                    i=0

    example:
           d(AAA, Dna) = 1 + 1 + 2 + 0 + 1 = 5
           ttaccttAAC   1
           gATAtctgtc   1
    Dna    ACGgcgttcg   2
           ccctAAAgag   0
           cgtcAGAggt   1

Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna) over all
k-mers Pattern, the same task that the Equivalent Motif Finding Problem is
trying to achieve. We call such a k-mer a median string for Dna.

Median String Problem
Find a median string.

Given: An integer k and a collection of strings Dna.

Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern.
(If multiple answers exist, you may return any one.)

Sample Dataset
3
AAATTGACGCAT
GACGACCACGTT
CGTCAGCGCCTG
GCTGAGCACCGG
AGTACGGGACAG

Sample Output
GAC

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> PREV: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference? (very important)
      -> HERE: "Median string problem" (BA2B) - not fast enough... find better algorithm
      ↓
    * GreedyMotifSearch
      -> NEXT: Profile-most probable k-mer problem (BA2C)

Info:
- Counting row-by-row returns exactly the same result as counting col-by-col.
                                                               (2)
                          T  C  G  G  G  G  g  T  T  T  t  t   3 row-by-row
                          c  C  G  G  t  G  A  c  T  T  a  C   4
                          a  C  G  G  G  G  A  T  T  T  t  C   2
                          T  t  G  G  G  G  A  c  T  T  t  t   4
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C   3
                          T  t  G  G  G  G  A  c  T  T  C  C   2
                          T  C  G  G  G  G  A  T  T  c  a  t   3
                          T  C  G  G  G  G  A  T  T  c  C  t   2
                          T  a  G  G  G  G  A  a  c  T  a  C   4
                          T  C  G  G  G  t  A  T  a  a  C  C + 3   Interesting!
    SCORE(Motifs)   (1)   3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30  result score is the same!
                          col-by-col

   (1) given motifs -> consensus
      ; Motifs -> CONSENSUS(Motifs)
      ; Motif Finding Problem
   (2) condidates of consensus -> best possible collection Motifs for each consensus string
      ; CONSENSUS(Motifs) -> Motifs
      ; Equivalent Motif Finding Problem

                      (1) first approach             (2) second approach
                      T C G G G G g T T T t t        T C G G G G g T T T t t   3
                      c C G G t G A c T T a C        c C G G t G A c T T a C   4
                      a C G G G G A T T T t C        a C G G G G A T T T t C   2
                      T t G G G G A c T T t t        T t G G G G A c T T t t   4
    Motifs            a a G G G G A c T T C C        a a G G G G A c T T C C   3
                      T t G G G G A c T T C C        T t G G G G A c T T C C   2
                      T C G G G G A T T c a t        T C G G G G A T T c a t   3
                      T C G G G G A T T c C t        T C G G G G A T T c C t   2
                      T a G G G G A a c T a C        T a G G G G A a c T a C   4
                      T C G G G t A T a a C C        T C G G G t A T a a C C + 3
    SCORE(Motifs)     3+4+0+0+1+1+1+5+2+3+6+4 = 30   3+4+0+0+1+1+1+5+2+3+6+4= 30

- comparing three problems
  ┌──────────────────────────────────────────────────────────────────┐
  │ Motif Finding Problem:                                           │ O(n^t * k * t)
  │   Given a collection of strings,                                 │
  │   find a set of k-mers, one from each string,                    │
  │   that minimizes the score of the resulting motif.               │
  │                                                                  │
  │ Input:                                                           │
  │   A collection of strings Dna and an integer k.                  │
  │ Output:                                                          │
  │   A collection "Motifs" of k-mers, one from each string in Dna,  │ Motifs
  │   minimizing SCORE(Motifs) among all possible choices of k-mers. │
  └──────────────────────────────────────────────────────────────────┘
  ┌────────────────────────────────────────────────────────────────────┐
  │ Equivalent Motif Finding Problem:                                  │
  │   Given a collection of strings,                                   │
  │   find a pattern and a collection of k-mers (one from each string) │
  │   that minimizes the distance between all possible patterns        │
  │   and all possible collections of k-mers.                          │
  │                                                                    │
  │ Input:                                                             │
  │   A collection of strings Dna and an integer k.                    │
  │ Output:                                                            │
  │   A k-mer "Pattern" and a collection of k-mers "Motifs",           │ Pattern, Motifs
  │   one from each string in Dna, minimizing d(Pattern, Motifs)       │
  │   among all possible choices of Pattern and Motifs.                │
  └────────────────────────────────────────────────────────────────────┘
  ┌──────────────────────────────────────────────────┐
  │ Median String Problem:                           │
  │   Find a median string.                          │
  │ Input:                                           │
  │   A collection of strings Dna and an integer k.  │
  │ Output:                                          │
  │   A k-mer Pattern minimizing d(Pattern, Dna)     │
  │   among all k-mers Pattern.                      │
  └──────────────────────────────────────────────────┘

- What is the meaning of above comparison?
  ; The col-by-col approach is exactly the same as "Motif finding problem", nothing else.
  ; But row-by-row approach is something different.
    The key observation for solving the Equivalent Motif Finding Problem is that,
    given Pattern, we don’t need to explore all possible collections (4^k) Motifs
    in order to minimize d(Pattern, Motifs).

    MOTIFS(AAA, Dna):
                         ttaccttAAC
                         gATAtctgtc
                     Dna ACGgcgttcg
                         ccctAAAgag
                         cgtcAGAggt

- What is the "median string"?
  ; a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern

- What are the meanings of "d(pattern, text)" and "d(pattern, Dna)"?
  d(pattern, one Dna string)       = x   (minimum hamming distance value)
  d(pattern, multiple Dna strings) = x + y + ... + z

Plan 1.
- pseudocode:
  ╔═══════════════════════════════════════════════════════╗
  ║  MEDIANSTRING(Dna, k)                                 ║ O(4^k * n * k * t)
  ║      distance <- infinite                             ║
  ║      for each k-mer Pattern from AA...AA to TT...TT   ║
  ║          if distance > d(Pattern, Dna)                ║
  ║              distance <- d(Pattern, Dna)              ║
  ║              Median <- Pattern                        ║
  ║      return Median                                    ║
  ╚═══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════
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

// BA01M: Number ro Pattern
fn number2pattern(mut number: usize, k: usize) -> String {
    let mut pattern: String = "".into();
    for i in 0..k {
        let temp = number / 4u32.pow((k - i - 1) as u32) as usize;
        pattern += match temp {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            _ => "-",
        };
        number -= temp * 4u32.pow((k - i - 1) as u32) as usize;
    }
    pattern
}

// d(pattern, text) function
fn min_distance_pattern_text(pattern: &str, text: &str) -> usize {
    let k: usize = pattern.len();
    let mut min_distance = usize::MAX;
    for i in 0..text.len() - k + 1 {
        let kmer = &text[i..i + k];
        let dist = hamming_distance(kmer, pattern);
        if dist < min_distance {
            min_distance = dist;
        }
    }
    return min_distance;
}

// d(pattern, texts) function
fn min_distance_pattern_texts(pattern: &str, texts: &[String]) -> usize {
    let mut min_distance = 0;
    for text in texts {
        min_distance += min_distance_pattern_text(pattern, text.as_str());
    }
    return min_distance;
}

// BA02B: median string
fn median_string(dnas: &[String], k: usize) -> HashSet<String> {
    let mut medians: HashSet<String> = HashSet::new();
    let mut min_distance: usize = usize::MAX;
    for idx in 0..4usize.pow(k as u32) - 1 {
        let pattern = number2pattern(idx, k);
        let dist = min_distance_pattern_texts(pattern.as_str(), dnas);
        if dist < min_distance {
            min_distance = dist;
            medians = HashSet::new();
            medians.insert(pattern.to_string());
        } else if dist == min_distance {
            medians.insert(pattern.to_string());
        }
    }
    return medians;
}

// BA02H: distance between pattern and strings
fn distance_bw_pattern_strings(pattern: &str, dnas: &Vec<String>) -> usize {
    let k: usize = pattern.len();
    let mut distance: usize = 0;
    for dna in dnas.iter() {
        let mut min_hdist = usize::MAX;
        for i in 0..dna.len() - k + 1 {
            let kmer = &dna[i..i + k];
            let hdistance = hamming_distance(pattern, kmer);
            if min_hdist > hdistance {
                min_hdist = hdistance;
            }
        }
        distance += min_hdist;
    }
    return distance;
}

// BA02H: median string
fn median_string_ba02h(dnas: &Vec<String>, k: usize) -> String {
    let mut min_distance = usize::MAX;
    let mut median = "".to_string();
    for i in 0..4usize.pow(k as u32) {
        let pattern = number2pattern(i, k);
        let distance = distance_bw_pattern_strings(&pattern, dnas);
        if min_distance > distance {
            min_distance = distance;
            median = pattern;
        }
    }
    return median;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02b.txt");
    let k = lines[0].parse::<usize>().unwrap();
    let dnas = &lines[1..];
    // dbg!(k, dnas);

    // as in textbook
    let start = Instant::now();
    for median_str in median_string(dnas, k) {
        println!("{}", median_str);
    }
    println!("{:?}", start.elapsed());

    // BA02H
    let start = Instant::now();
    println!("{}", median_string_ba02h(&dnas.to_vec(), k));
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
