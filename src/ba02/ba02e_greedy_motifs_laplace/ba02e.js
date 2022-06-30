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
  * Why do I need pseudo-count (Laplace’s Rule of Succession)?

                          A:  .2   .2  .0   .0   .0   .0   .9   .1   .1   .1   .3   .0
                          C:  .1   .6  .0   .0   .0   .0   .0   .4   .1   .2   .4   .6
                          G:  .0   .0   1    1   .9   .9   .1   .0   .0   .0   .0   .0
                          T:  .7   .2  .0   .0   .1   .1   .0   .5   .8   .7   .3   .4
  Pr(TCGTGGATTTCC|Profile) =  .7 * .6 * 1 * .0 * .9 * .9 * .9 * .5 * .8 * .7 * .4 * .6 = 0 (ZERO!!)

  * applying Laplace's Rule of Succession (Pseudo-count)

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

Plan 1.
  * modify 'PROFILE(Motifs)' function

═════════════════════════════════════════════════

References:
- Set a default parameter value for a JavaScript function
  https://stackoverflow.com/questions/894860/set-a-default-parameter-value-for-a-javascript-function
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

// Profile(Motifs), pseudocount
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

// BA02E: greedy motif search with pseudocount
function greedy_motif_search(dnas, k, t) {
    /**
     * ([str],int,int) -> [str]
     * returns a list of best motifs
     */
    // temporary best motifs
    let best_motifs = [];
    for (const dna of dnas) {
        best_motifs.push(dna.slice(0, k));
    }
    for (const kmer of all_kmers(dnas[0], k)) {
        let motifs = [kmer];
        for (let i = 1; i < t; i++) {
            const profile = get_profile_pseudocount(motifs);
            motifs.push([...profile_most_probable_kmer(dnas[i], k, profile)][0]);
        }
        if (get_score(motifs) < get_score(best_motifs)) {
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02e.txt').toString().split("\n");
    const [k, t] = lines[0].split(/\s/).map(Number);
    const dnas = lines.slice(1);

    const startTime = performance.now()
    const best_motifs = greedy_motif_search(dnas, k, t);
    for (const motif of best_motifs) {
        console.log(motif);
    }
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()