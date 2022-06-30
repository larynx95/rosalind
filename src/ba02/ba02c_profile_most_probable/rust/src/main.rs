/*
Rosalind: BA2C
Find a Profile-most Probable k-mer in a String

Given a profile matrix Profile,
we can evaluate the probability of every k-mer in a string Text
and find a Profile-most probable k-mer in Text,
i.e., a k-mer that was most likely to have been generated
by Profile among all k-mers in Text.
For example, ACGGGGATTACC is the Profile-most probable 12-mer in GGTACGGGGATTACCT.
Indeed, every other 12-mer in this string has probability 0.

In general, if there are multiple Profile-most probable k-mers in Text,
then we select the first such k-mer occurring in Text.

Given: A string Text, an integer k, and a 4 * k matrix Profile.

Return: A Profile-most probable k-mer in Text.
(If multiple answers exist, you may return any one.)

Sample Dataset
ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
5
0.2 0.2 0.3 0.2 0.3
0.4 0.3 0.1 0.5 0.1
0.3 0.3 0.5 0.2 0.4
0.1 0.2 0.1 0.1 0.2

Sample Output
CCGAG

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference? (very important)
      -> PREV: Median string problem (BA2B) - not fast enough... find another algorithm
      ↓
    * GreedyMotifSearch
      -> HERE: Profile-most probable k-mer problem (BA2C)
      -> NEXT: GreedyMotifSearch (BA2D)

Info:
    Profile-most Probable k-mer Problem
    Find a Profile-most probable k-mer in a string.
    (profile-most probable k-mer is not consensus.)

                          T  C  G  G  G  G  g  T  T  T  t  t
                          c  C  G  G  t  G  A  c  T  T  a  C
                          a  C  G  G  G  G  A  T  T  T  t  C
                          T  t  G  G  G  G  A  c  T  T  t  t
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C
                          T  t  G  G  G  G  A  c  T  T  C  C
                          T  C  G  G  G  G  A  T  T  c  a  t
                          T  C  G  G  G  G  A  T  T  c  C  t
                          T  a  G  G  G  G  A  a  c  T  a  C
                          T  C  G  G  G  t  A  T  a  a  C  C

    SCORE(Motifs)         3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30

    COUNT(Motifs)     A:  2  2  0  0  0  0  9  1  1  1  3  0
                      C:  1  6  0  0  0  0  0  4  1  2  4  6
                      G:  0  0 10 10  9  9  1  0  0  0  0  0
                      T:  7  2  0  0  1  1  0  5  8  7  3  4

    PROFILE(Motifs)   A: .2 .2  0  0  0  0 .9 .1 .1 .1 .3  0
                      C: .1 .6  0  0  0  0  0 .4 .1 .2 .4 .6
                      G:  0  0  1  1 .9 .9 .1  0  0  0  0  0
                      T: .7 .2  0  0 .1 .1  0 .5 .8 .7 .3 .4
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

                      A: .2 .2 .0 .0 .0 .0 .9 .1 .1 .1 .3 .0
                      C: .1 .6 .0 .0 .0 .0 .0 .4 .1 .2 .4 .6
                      G: .0 .0  1  1 .9 .9 .1 .0 .0 .0 .0 .0
                      T: .7 .2 .0 .0 .1 .1 .0 .5 .8 .7 .3 .4
    Pr(ACGGGGATTACC|Profile) = .2*.6*1*1*.9*.9*.9*.5*.8*.1*.4*.6 = 0.000839808
    (profile most probable kmer is not the same as concensus!!)

Plan 1.
- What is the "Greedy algorithm"?
  ; A greedy algorithm is an approach for solving a problem
    by selecting the best option available at the moment.
    It doesn't worry whether the current best result will bring the overall optimal result.
- Profile-most probable k-mer finding is a kind of greedy algorithm.

References:
- Create empty array with a length from a variable [duplicate]
  https://stackoverflow.com/questions/44847574/create-empty-array-with-a-length-from-a-variable
  Is it possible to have stack allocated arrays with the size determined at runtime in Rust?
  https://stackoverflow.com/questions/27859822/is-it-possible-to-have-stack-allocated-arrays-with-the-size-determined-at-runtim
- TODO: function argument: value or reference? Which is the better way?
- Creating two dimensional arrays in Rust
  https://stackoverflow.com/questions/13212212/creating-two-dimensional-arrays-in-rust
- How do I do insert/update of a Vec inside a HashMap
  https://users.rust-lang.org/t/how-do-i-do-insert-update-of-a-vec-inside-a-hashmap/17092
*/

// #!/usr/bin/env rust
use std::collections::HashMap;
use std::collections::HashSet;
use std::str::FromStr;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// get probability
fn get_probability(pattern: &str, profile: &Vec<Vec<f64>>) -> f64 {
    let mut pr: f64 = 1.0;
    let imap: HashMap<char, usize> = HashMap::from([('A', 0), ('C', 1), ('G', 2), ('T', 3)]);
    for (i, ch) in pattern.char_indices() {
        pr *= profile[imap[&ch]][i];
    }
    return pr;
}

// BA01C: profile most probable k-mer
fn profile_most_probable_kmers(text: &str, k: usize, profile: &Vec<Vec<f64>>) -> HashSet<String> {
    let mut answer: HashSet<String> = HashSet::new();
    let mut max_prob: f64 = f64::MIN;
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k];
        let prob = get_probability(pattern, profile);
        if prob > max_prob {
            max_prob = prob;
            answer = HashSet::new();
            answer.insert(pattern.to_string());
        } else if prob == max_prob {
            answer.insert(pattern.to_string());
        }
    }
    return answer;
}

// get probability, HashMap version
fn get_probability_hashmap(pattern: &str, profile: &HashMap<char, Vec<f64>>) -> f64 {
    let mut pr: f64 = 1.0;
    let imap: HashMap<char, usize> = HashMap::from([('A', 0), ('C', 1), ('G', 2), ('T', 3)]);
    for (i, ch) in pattern.char_indices() {
        pr *= profile[&ch][i];
    }
    return pr;
}

// BA01C: profile most probable k-mer, HashMap version
fn profile_most_probable_kmers_hashmap(
    text: &str,
    k: usize,
    profile: HashMap<char, Vec<f64>>,
) -> HashSet<String> {
    let mut answer: HashSet<String> = HashSet::new();
    let mut max_prob: f64 = f64::MIN;
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k];
        let prob = get_probability_hashmap(pattern, &profile);
        if prob > max_prob {
            max_prob = prob;
            answer = HashSet::new();
            answer.insert(pattern.to_string());
        } else if prob == max_prob {
            answer.insert(pattern.to_string());
        }
    }
    return answer;
}

// main
fn main() {
    // vector version
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02c.txt");
    let text = &lines[0];
    let k: usize = usize::from_str(&lines[1]).unwrap_or(0); // let k = lines[1].parse::<usize>().unwrap();
    let mut profile: Vec<Vec<f64>> = vec![];
    for line in &lines[2..] {
        // into_iter(), don't use it again
        let split = line
            .trim()
            .split(" ")
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        profile.push(split);
    }

    // hashmap version
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02c.txt");
    let text = &lines[0];
    let k: usize = usize::from_str(&lines[1]).unwrap_or(0);
    let mut temp_vec: Vec<Vec<f64>> = Vec::new();
    for line in &lines[2..] {
        // into_iter(), don't use it again
        let split = line
            .trim()
            .split(" ")
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        temp_vec.push(split);
    }
    let mut profile_hashmap: HashMap<char, Vec<f64>> = HashMap::new();
    for (i, nuc) in "ACGT".char_indices() {
        profile_hashmap.insert(nuc, temp_vec[i].clone());
    }

    // check time
    let start = Instant::now();
    for kmer in profile_most_probable_kmers(text, k, &profile) {
        println!("{} ", kmer);
    }
    println!("Execution time: {:?}", start.elapsed());

    // check time
    let start = Instant::now();
    for kmer in profile_most_probable_kmers_hashmap(text, k, profile_hashmap) {
        println!("{} ", kmer);
    }
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
