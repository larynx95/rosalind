/*
Rosalind: BA2H
Implement DistanceBetweenPatternAndStrings

The first potential issue with implementing MedianString from "Find a Median String" is
writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai),
the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}.
This task is achieved by the following pseudocode.

  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

Sample Dataset
AAA
TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

Sample Output
5

═════════════════════════════════════════════════

Pseudocode:
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ MEDIANSTRING(Dna, k)                                                   ║
  ║     distance 1                                                         ║
  ║     for i 0 to 4^k - 1                                                 ║
  ║         pattern <- NumberToPattern(i, k)                               ║
  ║         if distance > DistanceBetweenPaternAndStrings(Pattern, Dna)    ║
  ║             distance <- DistanceBetweenPatternAndStrings(Patern, Dna)  ║
  ║             Median <- Pattern                                          ║
  ║     return Median                                                      ║
  ╚════════════════════════════════════════════════════════════════════════╝

*/
// #!/usr/bin/env rust
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01G: hamming distance
fn hamming_distance(astr: &str, bstr: &str) -> usize {
    astr.char_indices()
        .filter(|&(i, ch)| ch != bstr.chars().nth(i).unwrap())
        .map(|(_, ch)| ch)
        .count()
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

// BA02H: distance between pattern and strings
fn distance_bw_pattern_strings(pattern: &str, dnas: &Vec<String>) -> usize {
    let k: usize = pattern.len();
    let mut distance: usize = 0;
    for dna in dnas.iter() {
        let mut min_hdist = usize::MAX;
        for i in 0..dna.len() - k + 1 {
            let kmer = &dna[i..i + k];
            let hdistance = hamming_distance(pattern, kmer);
            if min_hdist > hdistance {
                min_hdist = hdistance;
            }
        }
        distance += min_hdist;
    }
    return distance;
}

// BA02H: median string
fn median_string(dnas: &Vec<String>, k: usize) -> String {
    let mut min_distance = usize::MAX;
    let mut median = "".to_string();
    for i in 0..4usize.pow(k as u32) {
        let pattern = number2pattern(i, k);
        let distance = distance_bw_pattern_strings(&pattern, dnas);
        if min_distance > distance {
            min_distance = distance;
            median = pattern;
        }
    }
    return median;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02h.txt");
    let pattern = &lines[0];
    let dnas = lines[1]
        .trim()
        .split(" ")
        .map(|s| s.to_string())
        .collect::<Vec<String>>();
    //dbg!(pattern, split);

    // check time
    let start = Instant::now();
    println!("{}", distance_bw_pattern_strings(pattern, &dnas));
    println!("{:?}", start.elapsed());

    // check time
    let start = Instant::now();
    println!("{}", median_string(&dnas, 3));
    println!("{:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
