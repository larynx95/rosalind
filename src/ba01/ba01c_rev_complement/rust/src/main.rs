/*
Rosalind: BA1C
Find the Reverse Complement of a String

In DNA strings, symbols 'A' and 'T' are complements of each other,
as are 'C' and 'G'.
Given a nucleotide p, we denote its complementary nucleotide as p.
The reverse complement of a DNA string Pattern = p1…pn is the string rPattern = pn … p1
formed by taking the complement of each nucleotide in Pattern,
then reversing the resulting string.

For example, the reverse complement of Pattern = "GTCA" is rPattern = "TGAC".

Reverse Complement Problem
Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: rPattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT

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
      -> NEXT: Find Frequent Words with Mismatches and Reverse Complements (BA1)
        -> HERE: Find the Reverse Complement of a String (BA1C)

Plan 1.
- human DNA is double stranded helical structure, and has direction

  5-AAAACCCGGT-3
  3-  ...     -5     I should read this strand by 5 to 3 direction.

- steps:
  a. reverse given DNA strand
  b. get complement nucleotide for each character, that's all

═════════════════════════════════════════════════

References:
- Reversing a string in Rust
  https://stackoverflow.com/questions/27996430/reversing-a-string-in-rust
  let foo = "palimpsest";
  println!("{}", foo.chars().rev().collect::<String>());
- Original data is intact!
  Does 'map' function create a new iterator? Yes.
  fn map<B, F>(self, f: F) -> Map<Self, F>
  https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.map
  'map()' transforms one iterator into another, by means of its argument:
  something that implements 'FnMut'.
  It produces a new iterator which calls this closure on each element of the original iterator.
- Does 'rev' function create a new iterator, or mutate original one?
*/

use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
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

fn main() {
    // practice
    let str = "abcdef";
    // rev doesn't mutate original data
    println!("{}", str.chars().rev().collect::<String>());
    println!("{}", str);

    // map doesn't mutate original data
    println!("{}", str.chars().rev().map(|_| 'z').collect::<String>());
    println!("{}", str);

    // problem solving
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01c.txt");
    let pattern = &lines[0];
    dbg!(pattern);

    let start = Instant::now();
    println!("{}", reverse_complement(pattern));
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
