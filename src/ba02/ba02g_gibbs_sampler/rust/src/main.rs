/*
Rosalind: BA2G
Implement GibbsSampler

We have previously defined the notion of a Profile-most probable k-mer in a string.
We now define a Profile-randomly generated k-mer in a string Text.
For each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile),
resulting in n = |Text| - k + 1 probabilities (p1, …, pn).
These probabilities do not necessarily sum to 1,
but we can still form the random number generator Random(p1, …, pn) based on them.
GIBBSSAMPLER uses this random number generator
to select a Profile-randomly generated k-mer at each step:
if the die rolls the number i,
then we define the Profile-randomly generated k-mer as the i-th k-mer in Text.

Pseudocode:
  ╔═══════════════════════════════════════════════════════════════════════════════════╗
  ║ GIBBSSAMPLER(Dna, k, t, N)                                                        ║
  ║     randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna ║
  ║     BestMotifs <- Motifs                                                          ║
  ║     for j <- 1 to N                                                               ║
  ║         i <- Random(t)                                                            ║
  ║         Profile <- profile matrix constructed from all strings in Motifs          ║
  ║                    except for Motifi                                              ║
  ║         Motifi <- Profile-randomly generated k-mer in the i-th sequence           ║
  ║         if Score(Motifs) < Score(BestMotifs)                                      ║
  ║             BestMotifs <- Motifs                                                  ║
  ║     return BestMotifs                                                             ║
  ╚═══════════════════════════════════════════════════════════════════════════════════╝

Implement GibbsSampler
Given: Integers k, t, and N, followed by a collection of strings Dna.

Return: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts.
Remember to use pseudocounts!

Sample Dataset
8 5 100
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
      -> GreedyMotifSearch with pseudo-count (BA2E)
      ↓
    * Randomized Motif Search algorithm
      -> PREV: Randomized Motif Search (BA2F)
      -> HERE: Gibbs Sampling (BA2G)

Info:

    ttaccttAAC    tTACcttaac            ttaccttAAC    ttaccttAAC
    gATAtctgtc    gatATCtgtc            gATAtctgtc    gatatcTGTc
    ACGgcgttcg -> acggcgTTCg            ACGgcgttcg -> ACGgcgttcg
    ccctAAAgag    ccctaaAGAg            ccctAAAgag    ccctAAAgag
    cgtcAGAggt    CGTcagaggt            cgtcAGAggt    cgtcAGAggt

    RandomizedMotifSearch               GibbsSampler
    (may change all k-mers in one step) (change one-k-mer in one step)

  (1) first process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ttACCTtaac    taac     A: 2 1 1 1     A: 2/4 1/4 1/4 1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 1 1 1     C: 0   1/4 1/4 1/4
  Dna ccgGCGTtag  -> ---------- ->  /    -> G: 1 1 1 0  -> G: 1/4 1/4 1/4 0
      cactaACGAg     cactaACGAg    acta     T: 1 1 1 2     T: 1/4 1/4 1/4 2/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 2 2 2     A: 3/8 2/8 2/8 2/8
                      adding pseudocount -> C: 1 2 2 2  -> C: 1/8 2/8 2/8 2/8
                                            G: 2 2 2 1     G: 2/8 2/8 2/8 1/8
                                            T: 2 2 2 3     T: 2/8 2/8 2/8 3/8

    k-mers from deleted DNA string: ccgGCGTtag
    ccgG    cgGC    gGCG    GCGT    CGTt    GTta    Ttag
    0       0       0       1/128   0       1/256   0

    k-mers from deleted DNA string: ccgGCGTtag
    ccgG    cgGC    gGCG    *GCGT*  CGTt    GTta    Ttag
    4/8^4   8/8^4   8/8^4   24/8^4  12/8^4  16/8^4  8/8^4  total=80/8^4
    ─────   ─────   ─────   ──────  ──────  ──────  ─────
    80/8^4  80/8^4  80/8^4  80/8^4  80/8^4  80/8^4  80/8^4
=>  4/80    8/80    8/80    24/80   12/80   16/80   8/80   divided by total=80/8^4
=> RANDOM(4/80, 8/80, 8/80, 24/80, 12/80, 16/80, 8/80)     hypothetical seven-sided die

  (2) second process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ----------     /       A: 2 0 0 1     A: 2/4 0   0   1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 2 1 0     C: 0   2/4 1/4 0
  Dna ccgGCGTtag  -> ccgGCGTtag -> GCGT* -> G: 2 1 2 0  -> G: 2/4 1/4 2/4 0
      cactaACGAg     cactaACGAg    acta     T: 0 1 1 3     T: 0   1/4 1/4 3/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 1 1 2     A: 3/8 1/8 1/8 2/8
                      adding pseudocount -> C: 1 3 2 1  -> C: 1/8 3/8 2/8 1/8
                                            G: 3 2 3 1     G: 3/8 2/8 3/8 1/8
                                            T: 1 2 2 4     T: 1/8 2/8 2/8 4/8
    k-mers from deleted DNA string: ttACCTtaac
    ttAC    tACC    *ACCT*  CCTt    CTta    Ttaa    taac
    2/8^4   2/8^4   72/8^4  24/8^4  8/8^4   4/8^4   1/8^4

  (3) third process

      Dna                          Motifs  COUNT(Motifs)    PROFILE(Motifs)
      ttACCTtaac     ttACCTtaac    ACCT*    A: 2 0 0 1     A: 2/4 0   0   1/4
      gATGTctgtc     gATGTctgtc    GTct     C: 0 2 1 0     C: 0   2/4 1/4 0
  Dna ccgGCGTtag  -> ccgGCGTtag -> GCGT  -> G: 2 1 2 0  -> G: 2/4 1/4 2/4 0
      cactaACGAg     ----------     /       T: 0 1 1 3     T: 0   1/4 1/4 3/4
      cgtcagAGGT     cgtcagAGGT    AGGT
                                            A: 3 1 1 2     A: 3/8 1/8 1/8 2/8
                      adding pseudocount -> C: 1 3 2 1  -> C: 1/8 3/8 2/8 1/8
                                            G: 3 2 3 1     G: 3/8 2/8 3/8 1/8
                                            T: 1 2 2 4     T: 1/8 2/8 2/8 4/8
    k-mers from deleted DNA string: cactaACGAg
    cact    acta    ctaA    taAC    aACG    *ACGA*  CGAg
    15/8^4  9/8^4   2/8^4   1/8^4   9/8^4   27/8^4  2/8^4

═════════════════════════════════════════════════

References:
- enumerate
  https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.enumerate
- rand::distributions::weighted::WeightedIndex
  https://docs.rs/rand/latest/rand/distributions/weighted/struct.WeightedIndex.html
*/

use rand::distributions::WeightedIndex; // weighted rand sample
use rand::prelude::*;
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

// profile randomly generated k-mer
fn profile_randomly_generated_kmer(profile: &HashMap<char, Vec<f64>>, dna: &str) -> String {
    let k: usize = profile[&'A'].len(); // length of k-mer
    let n: usize = dna.len(); // length of dna string
                              // create a weighted vector with all k-mers in a dna string
    let mut rng = thread_rng();
    let mut wvec: Vec<(String, f64)> = Vec::new();
    for i in 0..n - k + 1 {
        let kmer = &dna[i..i + k];
        let prob = get_probability_hashmap(&kmer, profile);
        wvec.push((kmer.to_string(), prob));
    }
    // weighted random selection
    let dist = WeightedIndex::new(wvec.iter().map(|item| item.1)).unwrap();
    return wvec[dist.sample(&mut rng)].0.to_string();
}

// BA01G: Gibbs sampler
fn gibbs_sampler(dnas: &Vec<String>, k: usize, t: usize, n: usize) -> Vec<String> {
    // randomly select k-mers Motifs
    let mut rng = rand::thread_rng();
    let l: usize = dnas[0].len();
    let mut motifs: Vec<String> = Vec::new();
    for i in 0..t {
        let j: usize = rng.gen_range(0..l - k + 1);
        motifs.push(dnas[i][j..j + k].to_string());
    }
    let mut best_motifs: Vec<String> = motifs.clone();
    //
    for _ in 0..n {
        let rand_idx: usize = rng.gen_range(0..t);
        let rest = motifs
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx != rand_idx)
            .map(|(_, e)| e.to_string())
            .collect::<Vec<String>>();
        let profile = get_profile_pseudocount(&rest, 1);
        let motif_rand_idx = profile_randomly_generated_kmer(&profile, &dnas[rand_idx]);
        motifs[rand_idx] = motif_rand_idx;
        // compare score
        if get_score(&motifs) < get_score(&best_motifs) {
            best_motifs = motifs.clone();
        }
    }
    return best_motifs;
}

// repeat n times
fn repeat_n_times(dnas: &Vec<String>, k: usize, t: usize, n: usize, repeat: usize) -> Vec<String> {
    let mut best_score = usize::MAX;
    let mut best_motifs: Vec<String> = vec![];
    for _ in 0..repeat {
        let motifs = gibbs_sampler(dnas, k, t, n);
        let score = get_score(&motifs);
        if score < best_score {
            best_motifs = motifs;
            best_score = score;
        }
    }
    return best_motifs;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02g.txt");
    let split = lines[0]
        .trim()
        .split(" ")
        .map(|s| s.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    let k = split[0];
    let t = split[1];
    let n = split[2];
    let mut dnas: Vec<String> = Vec::new();
    for line in &lines[1..] {
        dnas.push(line.to_string());
    }
    // dbg!(dnas, k, t, n);

    // check time
    let start = Instant::now();
    for motif in repeat_n_times(&dnas, k, t, n, 20) {
        println!("{}", motif);
    }
    println!("Execution time: {:?}", start.elapsed()); // 227.2522659s, TODO: too slow!
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
