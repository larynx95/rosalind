/*
Rosalind: BA1H
Find All Approximate Occurrences of a Pattern in a String

We say that a k-mer Pattern appears as a substring of Text with at most d mismatches
if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern,
i.e., HammingDistance(Pattern, Pattern') <= d.
Our observation that a DnaA box may appear with slight variations
leads to the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return:
All starting positions where Pattern appears
as a substring of Text with at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3

Sample Output
6 7 26 27 78

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
      -> PREV: hamming distance (BA1G)
      -> HERE: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> frequency array
          -> NEXT: number to pattern (BA1M), pattern to number (BA1L)

Plan 1.
- pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  APPROXIMATEPATTERNCOUNT(Text, Pattern, d)        ║
  ║    count 0                                        ║
  ║    for i <- 0 to |Text| - |Pattern|               ║
  ║      Pattern' <- Text(i, |Pattern|)               ║
  ║      if HAMMINGDISTANCE(Pattern, Pattern') <= d   ║
  ║        count count + 1                            ║
  ║    return count                                   ║
  ╚═══════════════════════════════════════════════════╝
*/

use std::{
    fs::File,
    io::{prelude::*, BufReader},
    time::Instant,
};

// BA01G: hamming distance
fn hamming_distance(astr: &str, bstr: &str) -> usize {
    astr.char_indices()
        .filter(|&(i, ch)| ch != bstr.chars().nth(i).unwrap())
        .map(|(_, ch)| ch)
        .count()
}

// BA01H: All Approximate Occurrences (indices)
fn approximate_pattern_indices(text: &str, pattern: &str, d: usize) -> Vec<usize> {
    let k: usize = pattern.len();
    text.char_indices()
        .map(|(i, _)| i)
        .filter(|&i| i < text.len() - k + 1)
        .filter(|&i| hamming_distance(&text[i..(i + k)], pattern) <= d)
        .collect::<Vec<usize>>()
}

// BA01H: All Approximate Occurrences (count)
fn approximate_pattern_count(text: &str, pattern: &str, d: usize) -> usize {
    let k: usize = pattern.len();
    text.char_indices()
        .map(|(i, _)| i)
        .filter(|&i| i < text.len() - k + 1)
        .filter(|&i| hamming_distance(&text[i..(i + k)], pattern) <= d)
        .count()
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01h.txt");
    let pattern = &lines[0];
    let text = &lines[1];
    let d = lines[2].parse::<usize>().unwrap();

    let start = Instant::now();
    for index in approximate_pattern_indices(text, pattern, d) {
        print!("{index} ");
    }
    println!();
    println!("Execution time: {:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
