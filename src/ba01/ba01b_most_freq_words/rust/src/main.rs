/*
Rosalind: BA1B
Find the Most Frequent Words in a String

We say that Pattern is a most frequent k-mer in Text
if it maximizes Count(Text, Pattern) among all k-mers.
For example, "ACTAT" is a most frequent 5-mer
in "ACAACTATGCATCACTATCGGGAACTATCCT",
and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4

Sample Output
CATG GCAT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> HERE: find frequent words in a string (BA1B)
      -> NEXT: find all occurrence of a pattern in a string (BA1D)

Plan 1.
- writing down
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
- Pseudocode:
╔═════════════════════════════════════════════════╗
║ COMPUTINGFREQUENCIES(Text, k)                   ║
║   for i <- 0 to 4^k- 1                          ║
║     FREQUENCYARRAY(i) <- 0                      ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- Text(i, k)                       ║
║     j <- PATTERNTONUMBER(Pattern)               ║
║     FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1  ║
║   return FREQUENCYARRAY                         ║
╚═════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════╗
║ FREQUENTWORDS(Text, k)                          ║
║   FrequentPatterns an empty set                 ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- the k-mer Text(i, k)             ║
║     COUNT(i) <- PATTERNCOUNT(Text, Pattern)     ║
║   maxCount maximum value in array COUNT         ║
║   for i <- 0 to |Text| - k                      ║
║     if COUNT(i) = maxCount                      ║
║       add Text(i, k) to FrequentPatterns        ║
║   remove duplicates from FrequentPatterns       ║
║   return FrequentPatterns                       ║
╚═════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════╗
║ FASTERFREQUENTWORDS(Text , k)                     ║
║   FrequentPatterns <- an empty set                ║
║   FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k) ║
║   maxCount <- maximal value in FREQUENCYARRAY     ║
║   for i <- 0 to 4k - 1                            ║
║     if FREQUENCYARRAY(i) = maxCount               ║
║       Pattern <- NUMBERTOPATTERN(i, k)            ║
║       add Pattern to the set FrequentPatterns     ║
║   return FrequentPatterns                         ║
╚═══════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════╗
║ FINDINGFREQUENTWORDSBYSORTING(Text , k)              ║
║   FrequentPatterns <- an empty set                   ║
║   for i <- 0 to |Text| - k                           ║
║     Pattern <- Text(i, k)                            ║
║     INDEX(i) <- PATTERNTONUMBER(Pattern)             ║
║     COUNT(i) <- 1                                    ║
║   SORTEDINDEX <- SORT(INDEX)                         ║
║   for i <- 1 to |Text| - k                           ║
║     if SORTEDINDEX(i) = SORTEDINDEX(i - 1)           ║
║       COUNT(i) = COUNT(i - 1) + 1                    ║
║   maxCount <- maximum value in the array COUNT       ║
║   for i <- 0 to |Text| - k                           ║
║     if COUNT(i) = maxCount                           ║
║       Pattern <- NUMBERTOPATTERN(SORTEDINDEX(i), k)  ║
║       add Pattern to the set FrequentPatterns        ║
║   return FrequentPatterns                            ║
╚══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- Convert a String to int?
  https://stackoverflow.com/questions/27043268/convert-a-string-to-int
- pub fn parse<F>(&self) -> Result<F, <F as FromStr>::Err>
  https://doc.rust-lang.org/std/primitive.str.html#method.parse
- How do I get the key associated with the maximum value of a Rust HashMap?
  https://stackoverflow.com/questions/62525693/how-do-i-get-the-key-associated-with-the-maximum-value-of-a-rust-hashmap
- What is the difference between iter and into_iter?
  https://stackoverflow.com/questions/34733811/what-is-the-difference-between-iter-and-into-iter
*/

// #!/usr/bin/env rust

use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryInto;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03A: generate all k-mers
fn allkmers(text: &str, k: usize) -> Vec<&str> {
    let mut allkmers = Vec::new();
    for i in 0..text.len() - k + 1 {
        let kmer = &text[i..i + k];
        allkmers.push(kmer);
    }
    return allkmers;
}

// pattern count
fn pattern_count(text: &str, pattern: &str) -> usize {
    let mut cnt: usize = 0;
    let k = pattern.len() as usize;
    for i in 0..text.len() - k + 1 {
        let subs = &text[i..i + k];
        if subs == pattern {
            cnt += 1;
        }
    }
    return cnt;
}

// BA01B: frequent words
fn freq_words(text: &str, k: usize) -> HashSet<&str> {
    let mut fwords: HashSet<&str> = HashSet::new();
    let mut max_freq: usize = 0;
    for kmer in allkmers(text, k) {
        //
        let freq = pattern_count(text, kmer); //
        if freq > max_freq {
            max_freq = freq;
            fwords = HashSet::new();
            fwords.insert(kmer);
        } else if freq == max_freq {
            fwords.insert(kmer);
        }
    }
    return fwords;
}

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
    return num;
}

// BA01M: Number ro Pattern
fn number2pattern(mut number: u32, k: u32) -> String {
    let mut pattern: String = "".into();
    for i in 0..k {
        let temp = number / 4u32.pow((k - i - 1) as u32);
        pattern += match temp {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            _ => "-",
        };
        number -= temp * 4u32.pow((k - i - 1) as u32);
    }
    pattern
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
    return map;
}

// BA01B: faster frequent words
// TODO: borrow, ownership ... still confusing
fn freq_words_faster(text: &str, k: usize) -> HashSet<String> {
    let mut fwords: HashSet<String> = HashSet::new();
    let freq_map = compute_freq_hashmap(text, k);
    // let max_freq = freq_map.iter().max_by(|a, b| a.1.cmp(&b.1)).unwrap().1;  // borrowing freq_map
    let max_freq = freq_map.values().max().unwrap();
    for (key, value) in freq_map.iter() {
        if value == max_freq {
            let pattern = number2pattern(*key as u32, k as u32);
            fwords.insert(pattern);
        }
    }
    fwords
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01b.txt");
    let dna = &lines[0]; // :: &str
    let k = lines[1].parse::<usize>().unwrap();
    // dbg!(dna, k);

    // frequent words
    let start = Instant::now();
    for pattern in &freq_words(dna, k) {
        print!("{pattern} ");
    }
    println!();
    let elapsed = start.elapsed();
    println!("Execution time: {:?}", elapsed);

    // faster frequent words
    let start = Instant::now();
    for pattern in &freq_words_faster(dna, k) {
        print!("{pattern} ");
    }
    println!();
    let elapsed = start.elapsed();
    println!("Execution time: {:?}", elapsed);
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    // file automatically dropped when it is out-of scope.
    // no need to 'close' file
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
