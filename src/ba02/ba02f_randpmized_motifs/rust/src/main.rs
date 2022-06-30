/*
Rosalind: BA2F
Implement GreedyMotifSearch

We will now turn to randomized algorithms that flip coins and roll dice in order to search for motifs.
Making random algorithmic decisions may sound like a disastrous idea;
just imagine a chess game in which every move would be decided by rolling a die.
However, an 18th Century French mathematician and naturalist, Comte de Buffon,
first proved that randomized algorithms are useful by randomly dropping needles onto parallel strips of wood
and using the results of this experiment to accurately approximate the constant π.

Randomized algorithms may be nonintuitive because they lack the control of traditional algorithms.
Some randomized algorithms are Las Vegas algorithms,
which deliver solutions that are guaranteed to be exact,
despite the fact that they rely on making random decisions.
Yet most randomized algorithms are Monte Carlo algorithms.
These algorithms are not guaranteed to return exact solutions,
but they do quickly find approximate solutions.
Because of their speed, they can be run many times,
allowing us to choose the best approximation from thousands of runs.

A randomized algorithm for motif finding is given below.
  ╔══════════════════════════════════════════════════════════════════════════════════════╗
  ║  RANDOMIZEDMOTIFSEARCH(Dna, k, t)                                                    ║
  ║      randomly select k-mers Motifs = (Motif1, ... , Motift) in each string from Dna  ║
  ║      BestMotifs <- Motifs                                                            ║
  ║      while forever                                                                   ║
  ║          Profile <- Profile(Motifs)         # get Profile                            ║
  ║          Motifs <- Motifs(Profile, Dna)     # find best Motifs                       ║
  ║          if Score(Motifs) < Score(BestMotifs)                                        ║
  ║              BestMotifs <- Motifs                                                    ║
  ║          else                                                                        ║
  ║              return BestMotifs                                                       ║
  ╚══════════════════════════════════════════════════════════════════════════════════════╝

Implement RandomizedMotifSearch
Given: Positive integers k and t, followed by a collection of strings Dna.

Return: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times.
        Remember to use pseudocounts!

Sample Dataset
8 5
CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
TAGTACCGAGACCGAAAGAAGTATACAGGCGT
TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
AATCCACCAGCTCCACGTGCAATGTTGGCCTA

Sample Output
TCTCGGGG
CCAAGGTG
TACAGGCG
TTCAGGTG
TCCACGTG

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
      -> Median string problem (BA2B) - not fast enough... find another algorithm
      ↓
    * GreedyMotifSearch
      -> Profile-most probable k-mer problem (BA2C)
      -> GreedyMotifSearch (BA2D)
      -> PREV: GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> HERE: Randomized Motif Search (BA2F)
      -> NEXT: Gibbs Sampling (BA2G)

Chain of Functions:

  count ┌───────────────────────────────────────────────────────────── score ──┐
        └ get_profile ─── pr_score ─┬─ prof_most_probable_kmer ── get_motifs ──┼── random_motif_search
                         all_kmers ─┘                                          │
                                                               random_motifs ──┘

                     count
                 ┌─────┴────────────────┐
               get_profile            score
                 │                      │
  all_kmers    pr_score                 │
    └────────────┤                      │
         profile_most_probable_kmer     │    random_motifs
                 └──────────────────────┼──────┘
                             random_motif_search

Info:
- Review (what I have done until now)
  PROFILE(Motifs): get profile as a list of float list
  MOTIFS(Profile, Dnas): "profile-most probable k-mer"s from each DNA string

             A: 4/5   0    0   1/5        ttaccttaac
    profile  C:  0   3/5  1/5   0    DNA  gatgtctgtc
             G: 1/5  1/5  4/5   0         acggcgttag
             T:  0   1/5   0   4/5        ccctaacgag
                                          cgtcagaggt

    Motifs(Profile, Dna)   ttACCTtaac
                           gAGGTctgtc
                           acgGCGTtag
                           ccctaACGAg
                           cgtcagAGGT

- Monte-Carlo algorithm
    randomly chosen Motifs
    -> Motifs(Profile(Motifs), Dna)
    -> Profile(Motifs(Profile(Motifs), Dna))
    -> Motifs(Profile(Motifs(Profile(Motifs), Dna)))
    -> Profile(Motifs(Profile(Motifs(Profile(Motifs), Dna))))
    -> ... (a) get profile, (b) find best Motifs ... again and again

Plan 1.
- create a list of t random indices, then get random Motifs
  [idx1, idx, ... , idx_t]  -->  [DNA1[idx1:idx1+k], DNA2[idx2:idx2+k], ... , DNA_t[idx_t:idx_t + k]]

Plan 2.
- get all k-mers from each DNA strings, then randomly choose Motif from each collection of k-mers
  DNA1  --> k-mers --> Motif
  DNA2  --> k-mers --> Motif
  ...
  DNA_t --> k-mers --> Motif

═════════════════════════════════════════════════

References:
- Generate random numbers within a range
  https://rust-lang-nursery.github.io/rust-cookbook/algorithms/randomness.html
- Rust infinit loop
  loop
  https://doc.rust-lang.org/rust-by-example/flow_control/loop.html
*/

use rand::Rng; // 'rand' crate should be installed (add to Cargo.toml)
use std::collections::HashMap;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// 'Count' fx
fn get_count(motifs: &Vec<String>) -> HashMap<char, Vec<usize>> {
    let mut result: HashMap<char, Vec<usize>> =
        HashMap::from([('A', vec![]), ('C', vec![]), ('G', vec![]), ('T', vec![])]);
    // loop through rows
    for i in 0..motifs[0].len() {
        // loop through columns
        let mut temp: HashMap<char, usize> =
            HashMap::from([('A', 0), ('C', 0), ('G', 0), ('T', 0)]);
        for j in 0..motifs.len() {
            let nuc: char = motifs[j].chars().nth(i).unwrap();
            temp.insert(nuc, temp[&nuc] + 1);
        }
        // update result HashMap,
        for nuc in "ACGT".chars() {
            result.get_mut(&nuc).unwrap().push(temp[&nuc]) // get_mut
        }
    }
    return result;
}

// 'Score' fx
fn get_score(motifs: &Vec<String>) -> usize {
    let n: usize = motifs[0].len();
    let t: usize = motifs.len();
    let mut score: usize = 0;
    let count_map = get_count(motifs);
    for i in 0..n {
        let mut max_count: usize = 0;
        for ch in "ACGT".chars() {
            let freq: usize = count_map[&ch][i];
            if freq > max_count {
                max_count = freq;
            }
        }
        score += t - max_count;
    }
    return score;
}

// 'Profile' fx,
fn get_profile_pseudocount(motifs: &Vec<String>, pseudocount: usize) -> HashMap<char, Vec<f64>> {
    let t = motifs.len();
    let mut profile: HashMap<char, Vec<f64>> = HashMap::new();
    for (nuc, ivec) in get_count(motifs).iter() {
        let fvec: Vec<f64> = ivec
            .iter()
            .map(|&x| (x + pseudocount) as f64 / (t + t * pseudocount) as f64)
            .collect::<Vec<_>>();
        profile.insert(*nuc, fvec);
    }
    return profile;
}

// get probability, HashMap version
fn get_probability_hashmap(pattern: &str, profile: &HashMap<char, Vec<f64>>) -> f64 {
    let mut pr: f64 = 1.0;
    for (i, ch) in pattern.char_indices() {
        pr *= profile[&ch][i];
    }
    return pr;
}

// BA01C: profile most probable k-mer, HashMap version, return 'String'
fn profile_most_probable_kmers_hashmap(
    text: &str,
    k: usize,
    profile: &HashMap<char, Vec<f64>>,
) -> String {
    let mut answer: String = "".to_string();
    let mut max_prob: f64 = f64::MIN;
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k];
        let prob = get_probability_hashmap(pattern, &profile);
        if prob > max_prob {
            max_prob = prob;
            answer = pattern.to_string();
        }
    }
    return answer;
}

// get motifs
fn get_motifs(dnas: &Vec<String>, k: usize, profile: &HashMap<char, Vec<f64>>) -> Vec<String> {
    let mut motifs: Vec<String> = Vec::new();
    for dna in dnas {
        motifs.push(profile_most_probable_kmers_hashmap(dna, k, &profile));
    }
    return motifs;
}

// BA01F: randomized motif search
fn randomized_motif_search(dnas: &Vec<String>, k: usize, t: usize) -> Vec<String> {
    // get random motifs from each dna string
    let mut rng = rand::thread_rng();
    let n: usize = dnas[0].len();
    let mut motifs: Vec<String> = Vec::new();
    for i in 0..t {
        let j: usize = rng.gen_range(0..n - k + 1);
        motifs.push(dnas[i][j..j + k].to_string());
    }
    let mut best_motifs: Vec<String> = motifs.clone();
    // find best motifs
    loop {
        let profile: HashMap<char, Vec<f64>> = get_profile_pseudocount(&motifs, 1);
        motifs = get_motifs(&dnas, k, &profile);
        if get_score(&motifs) < get_score(&best_motifs) {
            best_motifs = motifs.clone();
        } else {
            return best_motifs;
        }
    }
}

// repeat n times
fn repeat_n_times(dnas: &Vec<String>, k: usize, t: usize, n: usize) -> Vec<String> {
    let mut min_score = usize::MAX;
    let mut best_motifs: Vec<String> = Vec::new();
    for _ in 0..n {
        let motifs: Vec<String> = randomized_motif_search(dnas, k, t);
        let score = get_score(&motifs);
        if score < min_score {
            min_score = score;
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02f.txt");
    let split = lines[0]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let t = split[1];
    let mut dnas: Vec<String> = Vec::new();
    for line in &lines[1..] {
        dnas.push(line.to_string());
    }
    // dbg!(dnas, k, t);

    // check time
    let start = Instant::now();
    for motif in repeat_n_times(&dnas, k, t, 1000) {
        println!("{}", motif);
    }
    println!("Execution time: {:?}", start.elapsed()); // 183.423 seconds!

    // TODO: Why is Rust so slow?
    // Rust      : 183.4235957         seconds
    // JavaScript:   2.697588200001046 seconds
    // Python    :  50.61851167678833  seconds
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
