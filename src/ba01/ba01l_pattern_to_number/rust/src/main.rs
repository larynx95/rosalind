/*
Rosalind: BA1L
Implement PatternToNumber

Implement PatternToNumber
Convert a DNA string to a number.

Given: A DNA string Pattern.
Return: PatternToNumber(Pattern).

Sample Dataset
AGT

Sample Output
11

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
      -> PREV: Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I) - brute-forced!
        -> NEXT: frequency array
          -> HERE: number to pattern (BA1M), pattern to number (BA1L)

Info.
  kmer    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
  index    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  x % 4    0  1  2  3  0  1  2  3  0  1  2  3  0  1  2  3
  x / 4    0  0  0  0  1  1  1  1  2  2  2  2  3  3  3  3
                A:0, C:1, G:2, T:3

    exp   1 0              2 1 0
    Nuc   T T : 3 3        G A T : 2 0 3
          │ └── 4^0*3      │ │ └── 4^0*2
          └──── 4^1*3      │ └──── 4^1*0
                 = 15      └────── 4^2*3
                                    = 50
Plan 1.
- pseudocode:
  ╔════════════════════════════════════════════════════════════════════╗
  ║  PATTERNTONUMBER(Pattern)                                          ║
  ║      if Pattern contains no symbols                                ║
  ║          return 0                                                  ║
  ║      symbol LASTSYMBOL(Pattern)                                    ║
  ║      Prefix PREFIX(Pattern)                                        ║
  ║      return 4 · PATTERNTONUMBER(Prefix) + SYMBOLTONUMBER(symbol)   ║
  ╚════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NUMBERTOPATTERN(index , k)                              ║
  ║      if k = 1                                            ║
  ║          return NUMBERTOSYMBOL(index)                    ║
  ║      prefixIndex QUOTIENT(index, 4)                      ║
  ║      r REMAINDER(index, 4)                               ║
  ║      symbol NUMBERTOSYMBOL(r)                            ║
  ║      PrefixPattern NUMBERTOPATTERN(prefixIndex, k  1)    ║
  ║      return concatenation of PrefixPattern with symbol   ║
  ╚══════════════════════════════════════════════════════════╝

  TODO: Is Recursion good choice in Rust programming language?
═════════════════════════════════════════════════

References:
- How to raise a number to a power?
  https://stackoverflow.com/questions/51208703/how-to-raise-a-number-to-a-power
- What is the default integer type in Rust? => i32
  https://stackoverflow.com/questions/55903243/what-is-the-default-integer-type-in-rust
- How to index a String in Rust
  https://stackoverflow.com/questions/24542115/how-to-index-a-string-in-rust

* TODO: Why does Rust's u64.pow expect a u32?
  https://stackoverflow.com/questions/57119562/why-does-rusts-u64-pow-expect-a-u32
  - pub const fn pow(self, exp: u32) -> u64
    https://doc.rust-lang.org/std/primitive.u64.html#method.pow
  - std::convert::TryInto
    https://doc.rust-lang.org/std/convert/trait.TryInto.html#tymethod.try_into

*/

use std::convert::TryInto;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01L: pattern to number
fn pattern2number(pattern: &str) -> usize {
    let mut num: usize = 0;
    let len = pattern.chars().count();
    // todo!();
    for i in 0..len {
        // let mult = 4u64.pow((len - i - 1) as u64);
        let mult = 4u64.pow(((len - i - 1) as u64).try_into().unwrap()) as usize;
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

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01l.txt");
    let pattern = &lines[0];
    // dbg!(pattern);

    let start = Instant::now();
    println!("{}", pattern2number(pattern));
    let elapsed = start.elapsed();
    println!("Execution time: {:?}", elapsed);
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
