/*
Rosalind: BA1D
Find All Occurrences of a Pattern in a String

In this problem, we ask a simple question:
how many times can one string occur as a substring of another?
Recall from “Find the Most Frequent Words in a String”
that different occurrences of a substring can overlap with each other.
For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
Find all occurrences of a pattern in a string.

Given:
Strings Pattern and Genome.

Return:
All starting positions in Genome where Pattern appears as a substring.
Use 0-based indexing.

Sample Dataset
ATAT
GATATATGCATATACTT

Sample Output
1 3 9

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
      -> HERE: find all occurrence of a pattern in a string (BA1D)
      -> NEXT: Clump finidng problem (BA1E)

Plan 1.
  * std::string::String::match_indices() fx can be used in overlapping pattern.
  * just use for-loop

═════════════════════════════════════════════════

Refrences:
- What is the correct way to check for string equality in JavaScript?
  https://stackoverflow.com/questions/3586775/what-is-the-correct-way-to-check-for-string-equality-in-javascript
*/

// #!/usr/bin/env rust
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01D: find all occurrences of a pattern in a text
fn find_all_indices_pattern(pattern: &str, genome: &str) -> Vec<usize> {
    let mut indices: Vec<usize> = vec![];
    let k = pattern.len();
    for i in 0..genome.len() - k + 1 {
        let subs = &genome[i..i + k];
        if subs == pattern {
            indices.push(i)
        }
    }
    return indices;
}

// BA01D: find all occurrences of a pattern in a text, filter version (method chaining)
fn find_all_indices_pattern_filter(pattern: &str, genome: &str) -> Vec<usize> {
    let k = pattern.len();
    genome
        .char_indices() // &str to iterator
        .map(|(i, _)| i) // extract indices
        .filter(|&i| i < genome.len() - k + 1) // get range of indices
        .filter(|&i| &genome[i..(i + k)] == pattern) // comparison  TODO: merge two filter steps
        .collect::<Vec<_>>() // iterator to Vector
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01d.txt");
    let pattern = &lines[0];
    let genome = &lines[1];

    let start = Instant::now();
    let indices = find_all_indices_pattern(pattern, genome);
    for index in indices {
        print!("{} ", index);
    }
    println!();
    println!("Execution time: {:?}", start.elapsed());

    let start = Instant::now();
    let indices = find_all_indices_pattern_filter(pattern, genome);
    for index in indices {
        print!("{} ", index);
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
