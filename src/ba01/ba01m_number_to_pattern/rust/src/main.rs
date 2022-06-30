/*
Rosalind: BA1M
Implement n2p

Convert an integer to its corresponding DNA string.

Given: Integers index and k.
Return: n2p(index, k).

Sample Dataset
45
4

Sample Output
AGTC

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

Plan 1.
- pseudocode:
  ╔════════════════════════════════════════════════════════════════════╗
  ║  PATTERNTONUMBER(Pattern)                                          ║
  ║      if Pattern contains no symbols                                ║
  ║          return 0                                                  ║
  ║      symbol LASTSYMBOL(Pattern)                                    ║
  ║      Prefix PREFIX(Pattern)                                        ║
  ║      return 4 · PATTERNTONUMBER(Prefix) + SYMBOLTONUMBER(symbol)   ║
  ║                                                                    ║
  ╚════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NUMBERTOPATTERN(index , k)                              ║
  ║      if k = 1                                            ║
  ║          return NUMBERTOSYMBOL(index)                    ║
  ║      prefixIndex <- QUOTIENT(index, 4)                   ║
  ║      r <- REMAINDER(index, 4)                            ║
  ║      symbol <- NUMBERTOSYMBOL(r)                         ║
  ║      PrefixPattern <- NUMBERTOPATTERN(prefixIndex, k-1)  ║
  ║      return concatenation of PrefixPattern with symbol   ║
  ║                                                          ║
  ╚══════════════════════════════════════════════════════════╝

References:
- Convert a String to int?
  https://stackoverflow.com/questions/27043268/convert-a-string-to-int
- How do you iterate over a string by character
  https://stackoverflow.com/questions/22118221/how-do-you-iterate-over-a-string-by-character
- Why is recursion not suggested in Rust?
  https://stackoverflow.com/questions/65948553/why-is-recursion-not-suggested-in-rust
*/

use std::{
    fs::File,
    io::{prelude::*, BufReader},
    time::Instant,
};

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

// BA01M: Number to Pattern
fn number2symbol(number: usize) -> String {
    match number {
        0 => "A".to_string(),
        1 => "C".to_string(),
        2 => "G".to_string(),
        3 => "T".to_string(),
        _ => "-".to_string(),
    }
}

fn number2pattern_rec(number: usize, k: usize) -> String {
    if k == 1 {
        return number2symbol(number);
    }
    let prefix_index = number / 4;
    let remains = number % 4;
    let prefix_pattern = number2pattern_rec(prefix_index, k - 1);
    return prefix_pattern + &number2symbol(remains);
}

fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01m.txt");
    let number = lines[0].parse::<usize>().unwrap();
    let k = lines[1].parse::<usize>().unwrap();
    // dbg!(number, k);

    let start = Instant::now();
    let pattern = number2pattern_rec(number, k);
    println!("{}", pattern);
    println!("Execution time: {:?}", start.elapsed());

    let start = Instant::now();
    let pattern = number2pattern(number, k);
    println!("{}", pattern);
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
