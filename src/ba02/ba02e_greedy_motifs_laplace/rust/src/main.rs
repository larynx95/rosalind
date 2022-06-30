/*
Rosalind: BA2E
Implement GreedyMotifSearch (pseudocount)

We encountered GreedyMotifSearch in “Implement GreedyMotifSearch”.
In this problem, we will power it up with pseudocounts.

Implement GreedyMotifSearch with Pseudocounts
Given: Integers k and t, followed by a collection of strings Dna.

Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t) with pseudocounts.
If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

Sample Dataset
3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG

Sample Output
TTC
ATC
TTC
ATC
TTC

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
      -> PREV: GreedyMotifSearch (BA2D)
      -> HERE: GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> NEXT: Randomized Motif Search (BA2F)

Info:
- Why do I need pseudo-count (Laplace’s Rule of Succession)?

                          A:  .2   .2  .0   .0   .0   .0   .9   .1   .1   .1   .3   .0
                          C:  .1   .6  .0   .0   .0   .0   .0   .4   .1   .2   .4   .6
                          G:  .0   .0   1    1   .9   .9   .1   .0   .0   .0   .0   .0
                          T:  .7   .2  .0   .0   .1   .1   .0   .5   .8   .7   .3   .4
  Pr(TCGTGGATTTCC|Profile) =  .7 * .6 * 1 * .0 * .9 * .9 * .9 * .5 * .8 * .7 * .4 * .6 = 0 (ZERO!!)

- applying Laplace’s Rule of Succession (Pseudo-count)
  Motifs            T    A    A    C
                    G    T    C    T
                    A    C    T    A
                    A    G    G    T

  COUNT(Motifs) A:  2    1    1    1       2+1   1+1   1+1   1+1
                C:  0    1    1    1  ->   0+1   1+1   1+1   1+1
                G:  1    1    1    0       1+1   1+1   1+1   0+1
                T:  1    1    1    2       1+1   1+1   1+1   2+1

  PROFILE(Motifs) 2/4  1/4  1/4  1/4       3/8   2/8   2/8   2/8
                    0  1/4  1/4  1/4  ->   1/8   2/8   2/8   2/8
                  1/4  1/4  1/4    0       2/8   2/8   2/8   1/8
                  1/4  1/4  1/4  2/4       2/8   2/8   2/8   3/8

═════════════════════════════════════════════════
References:
-
*/

use std::collections::HashMap;
use std::collections::HashSet;
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

// BA02E: greedy motif search with pseudocount
fn greedy_motif_search_pseudocount(dnas: &Vec<String>, k: usize, t: usize) -> Vec<String> {
    let mut best_motifs: Vec<String> = Vec::new();
    // get all first k-mers in each dna string
    for dna in dnas {
        best_motifs.push(dna[..k].to_string()); // == best_motifs.push(dna.chars().take(k).collect::<String>());
    }
    for i in 0..dnas[0].len() - k + 1 {
        // let mut motifs: Vec<String> = Vec::new();
        // motifs.push(dnas[0][i..i + k].to_string());
        let mut motifs: Vec<String> = Vec::from([dnas[0][i..i + k].to_string()]);
        for i in 1..t {
            let profile = get_profile_pseudocount(&motifs, 1);
            let motif_i = profile_most_probable_kmers_hashmap(dnas[i].as_str(), k, &profile);
            motifs.push(motif_i);
        }
        if get_score(&motifs) < get_score(&best_motifs) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02e.txt");
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
    for motif in greedy_motif_search_pseudocount(&dnas, k, t) {
        println!("{}", motif);
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
