/*
Rosalind: BA1I
Find the Most Frequent Words with Mismatches in a String

We defined a mismatch in "Compute the Hamming Distance Between Two Strings". We
now generalize "Find the Most Frequent Words in a String" to incorporate
mismatches as well.

Given strings Text and Pattern as well as an integer d, we define Count_d(Text,
Pattern) as the total number of occurrences of Pattern in Text with at most d
mismatches. For example, Count_1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because
AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA,
AAACA, and AAAGA. Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in Text is simply a string Pattern
maximizing Count_d(Text, Pattern) among all k-mers. Note that Pattern does not
need to actually appear as a substring of Text; for example, AAAAA is the most
frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA
does not appear exactly in this string. Keep this in mind while solving the
following problem.

Frequent Words with Mismatches Problem
Find the most frequent k-mers with mismatches in a string.

Given: A string Text as well as integers k and d.

Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
GATG ATGC ATGT

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
      -> HERE: Find the Most Frequent Words with Mismatches in a String (BA1I)
        -> frequency array
          -> number to pattern (BA1M), pattern to number (BA1L)
        -> PREV: neigbors (BA1N)

Plan 1.
- Pseudocode:
  ╔═══════════════════════════════════════════════════╗
  ║  APPROXIMATEPATTERNCOUNT(Text, Pattern, d)        ║
  ║    count 0                                        ║
  ║    for i 0 to |Text| - |Pattern|                  ║
  ║      Pattern <- Text(i, |Pattern|)                ║
  ║      if HAMMINGDISTANCE(Pattern, Pattern’) >= d   ║
  ║        count count + 1                            ║
  ║    return count                                   ║
  ╚═══════════════════════════════════════════════════╝

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

  ╔══════════════════════════════════════════════════════════════╗
  ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
  ║    Neighborhood <- set consisting of single string Pattern   ║
  ║    for j = 1 to d                                            ║
  ║      for each string Pattern' in Neighborhood                ║
  ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
  ║        remove duplicates from Neighborhood                   ║
  ║    return Neighborhood                                       ║
  ╚══════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  FREQUENTWORDSWITHMISMATCHES(Text, k, d)                              ║
  ║    FrequentPatterns an empty set                                      ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLOSE(i) <- 0                                                    ║
  ║      FREQUENCYARRAY <- 0                                              ║
  ║    for i <- 0 to |Text| - k                                           ║
  ║      Neighborhood <- NEIGHBORS(Text(i, k), d)                         ║
  ║      for each Pattern from Neighborhood                               ║
  ║        index <- p2n(Pattern)                                          ║
  ║        CLOSE(index) <- 1                                              ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLOSE(i) = 1                                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        FREQUENCYARRAY(i) <- APPROXIMATEPATTERNCOUNT(Text, Pattern, d) ║
  ║    maxCount maximal value in FREQUENCYARRAY                           ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if FREQUENCYARRAY(i) = maxCount                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════════════════╗
  ║  FINDINGFREQUENTWORDSWITHMISMATCHESBYSORTING(Text, k, d)                 ║
  ║    FrequentPatterns <- an empty set                                      ║
  ║    Neighborhoods <- an empty list                                        ║
  ║    for i <- 0 to |Text| - k                                              ║
  ║      add NEIGHBORS(Text(i, k), d) to Neighborhoods                       ║
  ║    form an array NEIGHBORHOODARRAY holding all strings in Neighborhoods  ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      Pattern <- NEIGHBORHOODARRAY(i)                                     ║
  ║      INDEX(i) <-  p2n(Pattern)                                           ║
  ║      COUNT(i) <- 1                                                       ║
  ║    SORTEDINDEX SORT(INDEX)                                               ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      if SORTEDINDEX(i) = SORTEDINDEX(i + 1)                              ║
  ║        COUNT(i + 1) <- COUNT(i) + 1                                      ║
  ║    maxCount maximum value in array COUNT                                 ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║     if COUNT(i) = maxCount                                               ║
  ║       Pattern <- n2p(SORTEDINDEX(i), k)                                  ║
  ║       add Pattern to FrequentPatterns                                    ║
  ║    return FrequentPatterns                                               ║
  ╚══════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════
References:
- How can I insert all values of one HashSet into another HashSet?
  https://stackoverflow.com/questions/56148586/how-can-i-insert-all-values-of-one-hashset-into-another-hashset
  a.extend(&b)
*/

// #!/usr/bin/env rust

use std::collections::HashSet;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
    time::Instant,
};

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

// BA01I: frequent words with mismatches
// TODO: This version is very slow...why? lots of repeating neighbors
fn freq_words_mismatch_slow(text: &str, k: usize, d: usize) -> HashSet<String> {
    let mut freq_patterns: HashSet<String> = HashSet::new();
    let mut max_count: usize = 0;
    for i in 0..text.len() - k + 1 {
        let neighborhood = recursive_neighbors(&text[i..i + k], d);
        // do multiple things all at once: count frequency, update max count, reset and insert to HashSet
        // But this part should be improved, because there's lots of unnecessary repeating neighbors!
        for nb in neighborhood {
            let count = approximate_pattern_count(text, &nb, d); // The number of neeighbors should be reduced!
            if count > max_count {
                max_count = count;
                freq_patterns = HashSet::from([nb]); // update max_count, reset HashSet, insert elem to HashSet
            } else if count == max_count {
                freq_patterns.insert(nb);
            }
        }
    }
    return freq_patterns;
}

// BA01J: frequent words with mismatches, improved
fn freq_words_mismatch(text: &str, k: usize, d: usize) -> HashSet<String> {
    let mut freq_patterns: HashSet<String> = HashSet::new();
    let mut max_count: usize = 0;
    let mut neighborhood: HashSet<String> = HashSet::new();
    // remove duplicates in neighborhood by HashSet
    for i in 0..text.len() - k + 1 {
        neighborhood.extend(recursive_neighbors(&text[i..i + k], d));
    }
    for nb in neighborhood {
        let count = approximate_pattern_count(text, &nb, d);
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
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01i.txt");
    let text = &lines[0];
    let split = lines[1]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let d = split[1];
    // dbg!(text, k, d);

    // slow version
    /*
    let start = Instant::now();
    for word in freq_words_mismatch_slow(text, k, d) {
        print!("{} ", word);
    }
    println!();
    println!("Execution time: {:?}", start.elapsed());
    */

    // frequent words mismatch reverse complement
    let start = Instant::now();
    for word in freq_words_mismatch(text, k, d) {
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
