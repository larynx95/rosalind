/*
Rosalind: BA1J
Find Frequent Words with Mismatches and Reverse Complements

We now extend “Find the Most Frequent Words with Mismatches in a String”
to find frequent words with both mismatches and reverse complements.
Recall that 'rPattern' refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
Find the most frequent k-mers (with mismatches and reverse complements)
in a DNA string.

Given: A DNA string Text as well as integers k and d.

Return:
All k-mers Pattern
maximizing the sum Count_d(Text, Pattern) + Count_d(Text, rPattern)
over all possible k-mers.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
ATGT ACAT

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
      -> Find the Most Frequent Words with Mismatches in a String (BA1I)
      -> HERE: Find Frequent Words with Mismatches and Reverse Complements (BA1J)
        -> Find the Reverse Complement of a String (BA1C)
      ↓
    * Solving the most frequent words with mismatch by improving algorithms
      -> Generate the Frequency Array of a String (BA1K)
        -> NEXT: pattern to number (BA1L), number to pattern (BA1M)

*/

// #!/usr/bin/env rust

use std::collections::HashSet;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
    time::Instant,
};

// BA01C: reverse complement
fn reverse_complement(dna: &str) -> String {
    dna.chars() // '&str' to interator
        .rev() //
        .map(|ch| match ch {
            // 'map' create a new iterator
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => 'X',
        })
        .collect::<String>()
}

// BA01G: hamming distance
fn hamming_distance(astr: &str, bstr: &str) -> usize {
    astr.char_indices()
        .filter(|&(i, ch)| ch != bstr.chars().nth(i).unwrap())
        .map(|(_, ch)| ch)
        .count()
}

// BA01N: recursive neighbors
fn recursive_neighbors(pattern: &str, d: usize) -> HashSet<String> {
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
    let suffix_neighbors = recursive_neighbors(&pattern[1..], d);
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

// BA01H: All Approximate Occurrences (count)
fn approximate_pattern_count(text: &str, pattern: &str, d: usize) -> usize {
    let k: usize = pattern.len();
    text.char_indices()
        .map(|(i, _)| i)
        .filter(|&i| i < text.len() - k + 1)
        .filter(|&i| hamming_distance(&text[i..(i + k)], pattern) <= d)
        .count()
}

// BA01J: frequent words with mismatches, reverse complement
fn freq_words_mismatch_RC(text: &str, k: usize, d: usize) -> HashSet<String> {
    let mut freq_patterns: HashSet<String> = HashSet::new();
    let mut max_count: usize = 0;
    let mut neighborhood: HashSet<String> = HashSet::new();
    // remove duplicates in neighborhood by HashSet
    for i in 0..text.len() - k + 1 {
        neighborhood.extend(recursive_neighbors(&text[i..i + k], d));
    }
    for nb in neighborhood {
        let count = approximate_pattern_count(text, &nb, d)
            + approximate_pattern_count(text, &reverse_complement(&nb), d);
        if count > max_count {
            max_count = count;
            freq_patterns = HashSet::from([nb]); // update max_count, reset HashSet, insert elem to HashSet
        } else if count == max_count {
            freq_patterns.insert(nb);
        }
    }
    return freq_patterns;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01j.txt");
    let text = &lines[0];
    let split = lines[1]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let d = split[1];
    // dbg!(text, k, d);

    // frequent words mismatch reverse complement
    let start = Instant::now();
    for word in freq_words_mismatch_RC(text, k, d) {
        print!("{} ", word);
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
