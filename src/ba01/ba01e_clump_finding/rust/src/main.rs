/*
Rosalind: BA1E
Find Patterns Forming Clumps in a String

Given integers L and t,
a string Pattern forms an (L, t)-clump inside a (larger) string Genome
if there is an interval of Genome of length L
in which Pattern appears at least t times.
For example, TGCA forms a (25,3)-clump in the following Genome:
gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem
Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset
CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4

Sample Output
CGACA GAAGA AATGT

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
      -> PREV: find all occurrence of a pattern in a string (BA1D)
      -> HERE: Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> NEXT: hamming distance (BA1G)

Plan 1.
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  CLUMPFINDING(Genome, k, t, L)                                        ║
  ║    FrequentPatterns <- an empty set                                   ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLUMP(i) <- 0                                                    ║
  ║    for i <- 0 to |Genome| - L                                         ║
  ║      Text <- the string of length L starting at position i in Genome  ║
  ║      FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)                  ║
  ║      for index <- 0 to 4k - 1                                         ║
  ║        if FREQUENCYARRAY(index) >= t                                  ║
  ║          CLUMP(index) <- 1                                            ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLUMP(i) = 1                                                  ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                               ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════╗
  ║  CLUMPFINDINGBetter(Genome, k, t, L)                      ║
  ║    FrequentPatterns <- an empty set                       ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      CLUMP(i) <- 0                                        ║
  ║    Text <- Genome(0, L)                                   ║
  ║    FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)        ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if FREQUENCYARRAY(i) >= t                            ║
  ║        CLUMP(i) <- 1                                      ║
  ║    for i <- 1 to |Genome| - L                             ║
  ║      FirstPattern <- Genome(i - 1, k)                     ║
  ║      index <- PATTERNTONUMBER(FirstPattern)               ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) - 1   ║
  ║      LastPattern <- Genome(i + L - k, k)                  ║
  ║      index <- PATTERNTONUMBER(LastPattern)                ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) + 1   ║
  ║      if FREQUENCYARRAY(index) >= t                        ║
  ║        CLUMP(index) <- 1                                  ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if CLUMP(i) = 1                                      ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                   ║
  ║        add Pattern to the set FrequentPatterns            ║
  ║    return FrequentPatterns                                ║
  ╚═══════════════════════════════════════════════════════════╝

1. find the generalized rule of frequency array

  "AAGCAAAGGTGGGC"  len=14
   vv-----------                                    first ... last
  "AAGCAAAGGTGGG"   len=13, starting from index 0   AA    ...  -
   "AGCAAAGGTGGGC"  len=13, starting from index 1   -     ...  GC
    -----------^^
    common part!

  (1) compute_freq("AAGCAAAGGTGGG", 2)
                    ^^  <--- minus 1
  (2) compute_freq( "AGCAAAGGTGGGC", 2)
                                ^^  <--- plus 1

      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
      3  0  2  0  1  0  0  0  0  1  3  1  0  0  1  0   --- (1)
     -2  0  2  0  1  0  0  0  0  2+ 3  1  0  0  1  0   --- (2)

2. If we know a frequency array of the first clump,
   we can get the frequency array of the whole genome.

    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT AAGCAAAGGTGGGC   fst  lst
    3  0  2  0  1  0  0  0  0  1  1  1  0  0  1  0  AAGCAAAGGTG
   -2  0  2  0  1  0  0  0  0  1  2+ 1  0  0  1  0   AGCAAAGGTGG     AA   GG
    2  0 -1  0  1  0  0  0  0  1  3+ 1  0  0  1  0    GCAAAGGTGGG    AG   GG
    2  0  1  0  1  0  0  0  0 -1+ 3  1  0  0  1  0     CAAAGGTGGGC   GC   GC
    check the frequency before subtraction!

═════════════════════════════════════════════════

References:
- How to iterate through a Hashmap, print the key/value and remove the value in Rust?
  https://stackoverflow.com/questions/45724517/how-to-iterate-through-a-hashmap-print-the-key-value-and-remove-the-value-in-ru
- How can I update a value in a mutable HashMap?
  https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap
  *my_map.get_mut("a").unwrap() += 10;
  *my_map.entry("a").or_insert(42) += 10;
- for and iterator
  https://doc.rust-lang.org/rust-by-example/flow_control/for.html#for-and-iterators
*/

// #!/usr/bin/env rust

use std::collections::HashMap;
use std::collections::HashSet; // HashSet
use std::convert::TryInto;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

/*
// BA03A: generate all k-mers
// to get L-length clump from genome, this fx can be useful
fn allkmers(text: &str, k: usize) -> Vec<&str> {
    // TODO: check the result type - Vec<&str> or Vec<String>, which is better?
    let mut allkmers = Vec::new();
    for i in 0..text.len() - k + 1 {
        let kmer = &text[i..i + k];
        allkmers.push(kmer);
    }
    return allkmers;
}
*/

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

// BA01E: clump finding
// for every L-clump, get frequency array, find patterns forming (L,t)-clump
fn clump_finding(genome: &str, k: usize, l: usize, t: usize) -> HashSet<String> {
    let mut freq_patterns: HashSet<String> = HashSet::new();
    for i in 0..genome.len() - l + 1 {
        let text = &genome[i..i + l];
        let freq_map = compute_freq_hashmap(text, k);
        for (key, value) in &freq_map {
            // &iterator
            if value >= &t {
                freq_patterns.insert(number2pattern(*key, k));
            }
        }
    }
    return freq_patterns;
}

// BA01E: clumpe finding, better
fn clump_finding_better(genome: &str, k: usize, l: usize, t: usize) -> HashSet<String> {
    let mut freq_patterns: HashSet<String> = HashSet::new();
    let mut freq_arr = compute_freq(&genome[0..l], k);
    // check the frequency of the first pattern in the first L-clump
    for idx in 0..freq_arr.len() {
        if freq_arr[idx] >= t {
            freq_patterns.insert(number2pattern(idx, k));
        }
    }
    // calculate the frequencies of the first and the last patterns, from left to right
    for i in 1..genome.len() - l + 1 {
        // updating the frequency of the first pattern
        let fst_pattern = &genome[i - 1..i - 1 + k];
        let idx = pattern2number(fst_pattern);
        freq_arr[idx] -= 1;
        // updating the frequency of the last pattern
        let lst_pattern = &genome[i + l - k..i + l];
        let idx = pattern2number(lst_pattern);
        freq_arr[idx] += 1;
        // and check the frequency of the last pattern
        if freq_arr[idx] >= t {
            freq_patterns.insert(number2pattern(idx, k));
        }
    }
    return freq_patterns;
}

// main
fn main() {
    // read a text file, assign values to variables
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01e.txt");
    let genome = &lines[0];
    let split = lines[1]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let l = split[1];
    let t = split[2];

    // check the result
    let start = Instant::now();
    for pattern in clump_finding(genome, k, l, t) {
        print!("{} ", pattern);
    }
    println!();
    println!("{:?}", start.elapsed());

    // check the result
    let start = Instant::now();
    for pattern in clump_finding_better(genome, k, l, t) {
        print!("{} ", pattern);
    }
    println!();
    println!("{:?}", start.elapsed());
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
