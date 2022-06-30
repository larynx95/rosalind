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
  * Review (what I have done until now)
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

  * Monte-Carlo algorithm
      randomly chosen Motifs
      -> Motifs(Profile(Motifs), Dna)
      -> Profile(Motifs(Profile(Motifs), Dna))
      -> Motifs(Profile(Motifs(Profile(Motifs), Dna)))
      -> Profile(Motifs(Profile(Motifs(Profile(Motifs), Dna))))
      -> ... (a) get profile, (b) find best Motifs ... again and again

Plan 1.
  * create a list of t random indices, then get random Motifs
    [idx1, idx, ... , idx_t]  -->  [DNA1[idx1:idx1+k], DNA2[idx2:idx2+k], ... , DNA_t[idx_t:idx_t + k]]

Plan 2.
  * get all k-mers from each DNA strings, then randomly choose Motif from each collection of k-mers
    DNA1  --> k-mers --> Motif
    DNA2  --> k-mers --> Motif
    ...
    DNA_t --> k-mers --> Motif

═════════════════════════════════════════════════

References:
- Getting a random value from a JavaScript array
  https://stackoverflow.com/questions/4550505/getting-a-random-value-from-a-javascript-array
  array[Math.floor(Math.random() * array.length)]
*/

// BA03A: all k-mers
function* all_kmers(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length - k + 1; i++) yield text.slice(i, i + k);
}

// Count(Motifs)
function get_count(motifs) {
    /**
     * [str] -> {char:[int]}
     * returns COUNT(Motifs)
     */
    // initialize dictionary
    const n = motifs[0].length;
    let dic_count = {};
    for (const nuc of ['A', 'C', 'G', 'T']) {
        dic_count[nuc] = [...Array(n).fill(0)];
    }
    // constructe dictionary
    for (const motif of motifs) {
        for (let i = 0; i < motif.length; i++) {
            switch (motif[i]) {
                case 'A':
                    dic_count['A'][i]++;
                    break;
                case 'C':
                    dic_count['C'][i]++;
                    break;
                case 'G':
                    dic_count['G'][i]++;
                    break;
                case 'T':
                    dic_count['T'][i]++;
                    break;
            }
        }
    }
    return dic_count;
}

// Score(Motifs)
function get_score(motifs) {
    /**
     * [str] -> int
     * returns SCORE(Motifs)
     */
    const n = motifs[0].length;
    const t = motifs.length;
    // create a 2d array of column-by-column COUNT(Motifs)
    let matrix_cbc = [...Array(n)].map(e => Array(4).fill(0));
    for (let i = 0; i < t; i++) {
        for (let j = 0; j < n; j++) {
            switch (motifs[i][j]) {
                case 'A':
                    matrix_cbc[j][0]++;
                    break;
                case 'C':
                    matrix_cbc[j][1]++;
                    break;
                case 'G':
                    matrix_cbc[j][2]++;
                    break;
                case 'T':
                    matrix_cbc[j][3]++;
                    break;
            }
        }
    }
    let score = 0;
    for (const col of matrix_cbc) {
        score += (t - Math.max.apply(Math, col));
    }
    return score
}

// Profile(Motifs) with pseudocount
function get_profile_pseudocount(motifs, pseudocount = 1) {
    /**
     * [str] -> {char:[float]}
     * returns a dictionary (object) of profile
     * embedding 'COUNT(Motifs)' function
     * >>> get_profile_pseudocount(["ACGTC", "ATTTT"])
     *     {A:[1,0,0,0,0],C:[0,0.5,0,0,0.5],G:[0,0,0.5,0,0],T:[0,0.5,0.5,1,0.5 ]}
     */
    const t = motifs.length;
    const n = motifs[0].length;
    let dic_profile = {};
    for (const nuc of ['A', 'C', 'G', 'T']) { dic_profile[nuc] = Array(n).fill(0); }
    for (const motif of motifs) {
        for (let i = 0; i < motif.length; i++) {
            switch (motif[i]) {
                case 'A':
                    dic_profile['A'][i]++;
                    break;
                case 'C':
                    dic_profile['C'][i]++;
                    break;
                case 'G':
                    dic_profile['G'][i]++;
                    break;
                case 'T':
                    dic_profile['T'][i]++;
                    break;
            }
        }
    }
    // profile by pseudocounting
    for (const nuc in dic_profile) {
        dic_profile[nuc] = dic_profile[nuc].map(num => (num + pseudocount) / (t + t * pseudocount));
    }
    return dic_profile
}

// Pr(pattern, profile)
function get_probability(pattern, profile) {
    /**
     * (str,{char:[float]}) -> float
     * returns a probability of a pattern
     */
    let pr = 1.0;
    for (let i = 0; i < pattern.length; i++) {
        switch (pattern[i]) {
            case 'A':
                pr *= profile['A'][i];
                break;
            case 'C':
                pr *= profile['C'][i];
                break;
            case 'G':
                pr *= profile['G'][i];
                break;
            case 'T':
                pr *= profile['T'][i];
                break;
        }
    }
    return pr;
}

// BA02C: profile most probable k-mer
function profile_most_probable_kmer(text, k, profile) {
    /**
     * (str,int,{char:[float]}) -> {str}
     * returns a set of profile most probable kmers
     */
    let result = new Set([]);
    let max_prob = -Infinity;
    for (let i = 0; i < text.length - k + 1; i++) {
        const kmer = text.slice(i, i + k);
        const prob = get_probability(kmer, profile);
        if (max_prob < prob) {
            max_prob = prob;
            result = new Set([kmer]);   // not 'Set(kmer)', Be careful.
        } else if (max_prob == prob) {
            result.add(kmer);
        }
    }
    return result;
}

// Motifs
function get_motifs(dnas, k, profile) {
    /**
     * ([str],int,{char:[float]}) -> [str]
     * returns a collection of Motifs
     */
    let motifs = [];
    for (const dna of dnas) {
        motifs.push([...profile_most_probable_kmer(dna, k, profile)][0]);
    }
    return motifs;
}

// BA02F: randomized motif search
function random_motif_search(dnas, k, t) {
    /**
     * ([str],int,int) -> [str]
     * returns a collection best motifs
     */
    let best_motifs = [];
    for (const dna of dnas) {
        const idx = Math.floor(Math.random() * (dna.length - k + 1));
        best_motifs.push(dna.slice(idx, idx + k));
    }
    while (true) {
        const profile = get_profile_pseudocount(best_motifs);
        const new_motifs = get_motifs(dnas, k, profile);
        if (get_score(new_motifs) < get_score(best_motifs)) {
            best_motifs = new_motifs;
        } else {
            return best_motifs;
        }
    }
}

// repeat n times
function repeat_manytimes(dnas, k, t, n) {
    let min_score = Infinity;
    let best_motifs = [];
    for (let i = 0; i < n; i++) {
        motifs = random_motif_search(dnas, k, t);
        const score = get_score(motifs);
        if (score < min_score) {
            min_score = score;
            best_motifs = motifs;
        }
    }
    return best_motifs
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02f.txt').toString().split("\n");
    const [k, t] = lines[0].split(/\s/).map(Number);
    const dnas = lines.slice(1);

    const startTime = performance.now()
    const motifs = repeat_manytimes(dnas, k, t, 1000);
    let answer = "";
    for (const motif of motifs) {
        console.log(motif);
    }
    console.log(`${performance.now() - startTime} milliseconds`)  // 144.6723958s
}

// execute main function
main()
