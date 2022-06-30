/*
Rosalind: BA1K
Generate the Frequency Array of a String

Given an integer k, we define the frequency array of a string Text as an array
of length 4k, where the i-th element of the array holds the number of times that
the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

kmer      AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
index      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
frequency  3   0   2   0   1   0   0   0   0   1   3   1   0   0   1   0

Computing a Frequency Array
Generate the frequency array of a DNA string.

Given: A DNA string Text and an integer k.

Return: The frequency array of k-mers in Text.

Sample Dataset
ACGCGGCTCTGAAA
2

Sample Output
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

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
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - by frequency array
        -> HERE: frequency array
          -> PREV: number to pattern (BA1M), pattern to number (BA1L)
        -> NEXT: neigbors (BA1N)

Plan 1.
- pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  COMPUTINGFREQUENCIES(Text, k)                    ║
  ║    for i <- 0 to 4^k- 1                           ║
  ║      FREQUENCYARRAY(i) <- 0                       ║
  ║    for i <- 0 to |Text| - k                       ║
  ║      Pattern <- Text(i, k)                        ║
  ║      j <- PATTERNTONUMBER(Pattern)                ║
  ║      FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1   ║
  ║    return FREQUENCYARRAY                          ║
  ╚═══════════════════════════════════════════════════╝

References:
- How do I return an array from a Rust function?
  https://stackoverflow.com/questions/59164456/how-do-i-return-an-array-from-a-rust-function
- Creating a vector of zeros for a specific size
  https://stackoverflow.com/questions/29530011/creating-a-vector-of-zeros-for-a-specific-size
- Can isize and usize be different in rust?
  https://stackoverflow.com/questions/55506647/can-isize-and-usize-be-different-in-rust
- Is there a method like JavaScript's substr in Rust?
  https://stackoverflow.com/questions/37157926/is-there-a-method-like-javascripts-substr-in-rust
  let s = "Hello, world!";  // :: &str
  let ss: String = s.chars().skip(7).take(5).collect();  // :: String
  let ss = &s[7..12];  // ::&str
- How can I update a value in a mutable HashMap?
  https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap
  - entry
    https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.entry
  - Enum std::collections::hash_map::Entry::or_insert()
    https://doc.rust-lang.org/std/collections/hash_map/enum.Entry.html#method.or_insert
*/

use std::collections::HashMap;
use std::convert::TryInto;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
    time::Instant,
}; // dictionary in Python

// BA01L: pattern to number
fn pattern2number(pattern: &str) -> u64 {
    let mut num: u64 = 0;
    let len = pattern.chars().count();
    // todo!();
    for i in 0..len {
        // let mult = 4u64.pow((len - i - 1) as u64);
        let mult = 4u64.pow(((len - i - 1) as u64).try_into().unwrap()); // TODO: u64.pow(u32)
        let num_char = match pattern.as_bytes()[i] as char {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => 0,
        };
        num += mult * num_char;
    }
    num
}

// BA01K: compute frequency
fn compute_freq(text: &str, k: usize) -> Vec<usize> {
    let mut freq_vec = vec![0; 4u32.pow(k as u32) as usize];
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k]; // :: &str
        let number = pattern2number(pattern);
        freq_vec[number as usize] += 1;
    }
    return freq_vec;
}

// BA01K: coumpute frequency, HashMap version
fn compute_freq_hashmap(text: &str, k: usize) -> HashMap<usize, usize> {
    let mut map = HashMap::new(); // TODO: How did Rust know the type of map?
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k]; // :: &str
        let idx = pattern2number(pattern) as usize;
        *map.entry(idx).or_insert(0) += 1; // if not exists, insert 0, then +1
                                           // if exists, then +1
    }
    map
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01k.txt");
    let dna = &lines[0];
    let k = lines[1].parse::<usize>().unwrap();
    // dbg!(dna, k);

    // vector version
    let start = Instant::now();
    for i in compute_freq(dna, k) {
        print!("{i} ");
    }
    println!();
    let elapsed = start.elapsed();
    println!("Execution time: {:?}", elapsed); // 53.7 us

    // HashMap version
    let start = Instant::now();
    let map = compute_freq_hashmap(dna, k);
    for i in 0..4u32.pow(k as u32) as usize {
        if map.contains_key(&i) {
            print!("{} ", map[&i]);
        } else {
            print!("0 ");
        }
    }
    println!();
    let elapsed = start.elapsed();
    println!("Execution time: {:?}", elapsed); // 72.6 us
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
