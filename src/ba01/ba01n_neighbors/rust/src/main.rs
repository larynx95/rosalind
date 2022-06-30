/*
Rosalind: BA1N
Generate the d-Neighborhood of a String

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming
distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String
Find all the neighbors of a pattern.

Given: A DNA string Pattern and an integer d.
Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
ACG
1

Sample Output
CCG
TCG
GCG
AAG
ATG
AGG
ACA
ACC
ACT
ACG

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
      -> NEXT: Find the Most Frequent Words with Mismatches in a String (BA1I)
        -> PREV: frequency array
          -> number to pattern (BA1M), pattern to number (BA1L)
        -> HERE: neigbors (BA1N)

Plan 1.
  * immediate neighbors
    ACG, d=1
    A --> CCG, GCG, TCG
    C --> AAG, AGG, ATG
    G --> ACA, ACC, ACT
    ╔═══════════════════════════════════════════════════════════════════════╗
    ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
    ║    Neighborhood <- the set consisting of the single string Pattern    ║
    ║    for i = 1 to |Pattern|                                             ║
    ║      symbol <- i-th nucleotide of Pattern                             ║
    ║      for each nucleotide x different from symbol                      ║
    ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
    ║        add Neighbor to Neighborhood                                   ║
    ║    return Neighborhood                                                ║
    ╚═══════════════════════════════════════════════════════════════════════╝

  * iterative neighbors
    ╔══════════════════════════════════════════════════════════════╗
    ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
    ║    Neighborhood <- set consisting of single string Pattern   ║
    ║    for j = 1 to d                                            ║
    ║      for each string Pattern' in Neighborhood                ║
    ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
    ║        remove duplicates from Neighborhood                   ║
    ║    return Neighborhood                                       ║
    ╚══════════════════════════════════════════════════════════════╝

  * neighbors - recursive function
    (1) base case
      - if d == 0          => {pattern}
      - if |pattern| == 0  => {}
      - if |pattern| == 1  => {'A', 'C', 'G', 'T'}
    (2) recursive case
      - if neighbors(pattern[1:]) < d   => can prepend any nucleotide
      - if neighbors(pattern[1:]) >= d  => can change the first nucleotide
    ╔══════════════════════════════════════════════════════════╗
    ║  NEIGHBORS(Pattern, d)                                   ║
    ║    if d = 0                                              ║
    ║      return {Pattern}                                    ║
    ║    if |Pattern| = 1                                      ║
    ║      return {A, C, G, T}                                 ║
    ║    Neighborhood <- an empty set                          ║
    ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
    ║    for each string Text from SuffixNeighbors             ║
    ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
    ║        for each nucleotide x                             ║
    ║          add x + Text to Neighborhood                    ║
    ║      else                                                ║
    ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
    ║    return Neighborhood                                   ║
    ╚══════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- How do I change characters at a specific index within a string in rust?
  https://stackoverflow.com/questions/66661118/how-do-i-change-characters-at-a-specific-index-within-a-string-in-rust
- replace_range
  https://doc.rust-lang.org/std/string/struct.String.html#method.replace_range
- str::char_indices
  https://doc.rust-lang.org/std/primitive.str.html#method.char_indices
- How can I insert all values of one HashSet into another HashSet?
  https://stackoverflow.com/questions/56148586/how-can-i-insert-all-values-of-one-hashset-into-another-hashset
  https://doc.rust-lang.org/std/iter/trait.Extend.html#tymethod.extend
- How can I do a mutable borrow in a for loop?
  https://stackoverflow.com/questions/39622783/how-can-i-do-a-mutable-borrow-in-a-for-loop
- Why does HashMap have iter_mut() but HashSet doesn't?
  https://stackoverflow.com/questions/35970238/why-does-hashmap-have-iter-mut-but-hashset-doesnt

*/

// #!/usr/bin/env rust

use std::collections::HashSet;
use std::{
    fs::File,
    io::{prelude::*, BufReader, Write},
    time::Instant,
};

// BA01N: immediate neighbors
fn immediate_neighbors(pattern: &str) -> HashSet<String> {
    let mut neighborhood: HashSet<String> = HashSet::new();
    neighborhood.insert(pattern.to_string());
    for (i, ch) in pattern.char_indices() {
        // char_indices() create iterator (tuple)
        for nuc in "ACGT".chars() {
            // chars() create iterator
            if ch != nuc {
                let mut neighbor = pattern.clone().to_string(); // replace_range() is destructive
                neighbor.replace_range(i..i + 1, &nuc.to_string()); // need &str, not String
                neighborhood.insert(neighbor);
            }
        }
    }
    return neighborhood;
}

// BA01N: iterative neighbors
// TODO: Notice this for-loop!
// Each time a for-loop is repeated, the set changes!
// Is this possible?
// Without clone(), this for-loop does not work at all.
fn iterative_neighbors(pattern: &str, d: usize) -> HashSet<String> {
    let mut neighborhood: HashSet<String> = HashSet::from([pattern.to_string()]);
    for _ in 0..d {
        for neighbor in neighborhood.clone() {
            neighborhood.extend(immediate_neighbors(&neighbor));
        }
    }
    return neighborhood;
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

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01n.txt");
    let pattern = &lines[0];
    let k = lines[1].parse::<usize>().unwrap();

    // immediate neighbors
    let start = Instant::now();
    let neighborhood = immediate_neighbors(pattern);
    for nb in &neighborhood {
        // by reference
        println!("{}", nb);
    }
    println!("{}", neighborhood.len()); // by value
    println!("Execution time (immediate): {:?}", start.elapsed());

    // iterative neighbors
    let start = Instant::now();
    let neighborhood = iterative_neighbors(pattern, k);
    for nb in &neighborhood {
        println!("{}", nb);
    }
    println!("{}", neighborhood.len());
    println!("Execution time (itrative): {:?}", start.elapsed());

    // recursive neighbors
    let start = Instant::now();
    let neighborhood = recursive_neighbors(pattern, k);
    for nb in &neighborhood {
        println!("{}", nb);
    }
    println!("{}", neighborhood.len());
    println!("Execution time (recursive): {:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}

// helper function: write to text file
fn answer_to_file(filename: &str, result: HashSet<String>) {
    // write to text file
    let mut w = File::create(filename).unwrap();
    for elem in result {
        writeln!(&mut w, "{}", elem).unwrap();
    }
}
