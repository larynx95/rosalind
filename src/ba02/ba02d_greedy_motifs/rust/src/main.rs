/*
Rosalind: BA2D
Implement GreedyMotifSearch

  ╔═════════════════════════════════════════════════════════════════════════════════╗
  ║  GREEDYMOTIFSEARCH(Dna, k, t)                                                   ║
  ║      BestMotifs <- motif matrix formed by first k-mers in each string from Dna  ║
  ║      for each k-mer Motif in the first string from Dna                          ║
  ║          Motif1 <- Motif                                                        ║
  ║          for i = 2 to t                                                         ║
  ║              form Profile from motifs Motif1, ..., Motif_i-1                    ║
  ║              Motifi <- Profile-most probable k-mer in the i-th string in Dna    ║
  ║          Motifs <- (Motif1, ..., Motift)                                        ║
  ║          if SCORE(Motifs) < SCORE(BestMotifs)                                   ║
  ║              BestMotifs <- Motifs                                               ║
  ║      return BestMotifs                                                          ║
  ╚═════════════════════════════════════════════════════════════════════════════════╝

Implement GreedyMotifSearch
Given: Integers k and t, followed by a collection of strings Dna.

Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t).
If at any step you find more than one Profile-most probable k-mer in a given string,
use the one occurring first.

Sample Dataset
3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG

Sample Output
CAG
CAG
CAA
CAA
CAA

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
      -> PREV: Profile-most probable k-mer problem (BA2C)
      -> HERE: GreedyMotifSearch (BA2D)
      -> NEXT: GreedyMotifSearch with pseudo-count (BA2E)

Info.
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
-
- implementing GreedyMotifSearch using "Profile-most probable k-mer" function in BA2C
- Pseudocode:
  ╔═════════════════════════════════════════════════════════════════════════════════╗
  ║  GREEDYMOTIFSEARCH(Dna, k, t)                                                   ║
  ║      BestMotifs <- motif matrix formed by first k-mers in each string from Dna  ║
  ║      for each k-mer Motif in the first string from Dna                          ║
  ║          Motif1 <- Motif                                                        ║
  ║          for i = 2 to t                                                         ║
  ║              form Profile from motifs Motif1, ..., Motif_i-1                    ║
  ║              Motifi <- Profile-most probable k-mer in the i-th string in Dna    ║
  ║          Motifs <- (Motif1, ..., Motift)                                        ║
  ║          if SCORE(Motifs) < SCORE(BestMotifs)                                   ║
  ║              BestMotifs <- Motifs                                               ║
  ║      return BestMotifs                                                          ║
  ╚═════════════════════════════════════════════════════════════════════════════════╝

- practice with sample dataset
  not easy, divide algorithm into several steps
  (1) Let's make two assumptions.
      One for "best Motifs", the other for best "Motif"
      Motifs = {Motif1, Motif2, ... , Motif_t}

      Assumption 1. The best Motifs == first k-mers in each DNA strings (len(Motifs) == t)
      Assumption 2. The best Motif in a DNA string == first k-mer

                              ↓first                              ↓last
      GGCGTTCAGGCA --> kmers: GGC GCG CGT GTT TTC TCA CAG AGG GGC GCA
      AAGAATCAGTCA            one of these k-mers will be a Motif in the first DNA string
      CAAGGAGTTCGC            Let's suppose the Motif is the first k-mer 'GGC' (by Assumption 1.),
      CACGTCAATCAC            and best Motifs is ['GGC','AAG','CAA','CAC','CAA'].
      CAATAATATTCG

  (2) Put the best Motifs (Motifs_0: first t k-mers (t * Motif_0)) aside for a while.

  (3) Find the best Motifs (Motifs_i) for each i-th k-mer (Motif_i) in the first DNA string.
      And get the NEW best Motifs by comparing SCORE(Motifs_0) with SCORE(Motifs_i).
      Do this process again and again till the last k-mer in the first DNA string.
      I can get the real BEST Motifs by this "len(DNA[0]) - k + 1" iteration.

- TODO: GreedyMotifSearch fx in BA02D is not always correct.
  Because 'Profile-most probable k-mer' fx can have more than one result.
  To pass this BA02D, take the first k-mer in 'Profile-most probable k-mer',
  not select one from the result HashSet.

═════════════════════════════════════════════════

References:
* Greedy Algorithm
  https://www.programiz.com/dsa/greedy-algorithm#:~:text=A%20greedy%20algorithm%20is%20an,if%20the%20choice%20is%20wrong.
* modifying vector value in HashMap
  - How can I update a value in a mutable HashMap?
    https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap
  - HashMap or_insert()
    https://doc.rust-lang.org/std/collections/hash_map/enum.Entry.html#method.or_insert
  - How to build a HashMap of Vectors in Rust? (uncomfortable)
    https://stackoverflow.com/questions/26169216/how-to-build-a-hashmap-of-vectors-in-rust
- What is the difference between iter and into_iter?
  https://stackoverflow.com/questions/34733811/what-is-the-difference-between-iter-and-into-iter
- Using map with Vectors
  https://stackoverflow.com/questions/30026893/using-map-with-vectors
- If only one element in a hashset, how can I get it out?
  https://stackoverflow.com/questions/23595749/if-only-one-element-in-a-hashset-how-can-i-get-it-out
  set.iter().next().unwrap();
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
        // update result HashMap, TODO: This part is very difficult.
        for nuc in "ACGT".chars() {
            result.get_mut(&nuc).unwrap().push(temp[&nuc])
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
fn get_profile(motifs: &Vec<String>) -> HashMap<char, Vec<f64>> {
    let t = motifs.len();
    let mut profile: HashMap<char, Vec<f64>> = HashMap::new();
    for (nuc, ivec) in get_count(motifs).iter() {
        let fvec: Vec<f64> = ivec
            .iter()
            .map(|&x| x as f64 / t as f64)
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

// BA01C: profile most probable k-mer, HashMap version
// There may be more than one result.
// But in BA02D, take the first one.
fn profile_most_probable_kmers_hashmap_ret_hashset(
    text: &str,
    k: usize,
    profile: HashMap<char, Vec<f64>>,
) -> HashSet<String> {
    let mut answer: HashSet<String> = HashSet::new();
    let mut max_prob: f64 = f64::MIN;
    for i in 0..text.len() - k + 1 {
        let pattern = &text[i..i + k];
        let prob = get_probability_hashmap(pattern, &profile);
        // TODO: In BA02D, just take the first one, not all k-mers.
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

// BA01C: profile most probable k-mer, HashMap version, return 'String'
// In BA02D, I modified the return type from 'HashSet<String>' to 'String'.
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

// BA02D: greedy motif search
fn greedy_motif_search(dnas: &Vec<String>, k: usize, t: usize) -> Vec<String> {
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
            let profile = get_profile(&motifs);
            let motif_i = profile_most_probable_kmers_hashmap(dnas[i].as_str(), k, &profile);
            motifs.push(motif_i);
        }
        if get_score(&motifs) < get_score(&best_motifs) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba02d.txt");
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
    // dbg!(&dnas, k, t);

    // check time
    let start = Instant::now();
    for motif in greedy_motif_search(&dnas, k, t) {
        println!("{}", motif);
    }
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
