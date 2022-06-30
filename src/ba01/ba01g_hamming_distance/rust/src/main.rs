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

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01g.txt");
    let astr = &lines[0];
    let bstr = &lines[1];

    let start = Instant::now();
    println!("{}", hamming_distance(astr, bstr));
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
