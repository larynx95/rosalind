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

use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03A: get all k-mers
fn all_kmers(text: &str, k: usize) -> Vec<String> {
    let length = text.len();
    let mut kmers: Vec<String> = Vec::new();
    for (i, _) in text.char_indices().take(length - k + 1) {
        kmers.push(text[i..(i + k)].to_string());
    }
    return kmers;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03a.txt");
    let k: usize = lines[0].parse::<usize>().unwrap();
    let text = &lines[1];

    // check time
    let start = Instant::now();
    for pat in &all_kmers(text, k) {
        println!("{}", pat);
    }
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
