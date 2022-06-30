/*
Rosalind: BA2D
Implement GreedyMotifSearch

╔════════════════════════════════════════════════════════════════════════════════╗
║ GREEDYMOTIFSEARCH(Dna, k, t)                                                   ║
║     BestMotifs <- motif matrix formed by first k-mers in each string from Dna  ║
║     for each k-mer Motif in the first string from Dna                          ║
║         Motif_1 <- Motif                                                       ║
║         for i = 2 to t                                                         ║
║             form Profile from motifs Motif_1, ..., Motif_i - 1                 ║
║             Motif_i <- Profile-most probable k-mer in the i-th string in Dna   ║
║         Motifs <- (Motif_1, ..., Motif_t)                                      ║
║         if Score(Motifs) < Score(BestMotifs)                                   ║
║             BestMotifs <- Motifs                                               ║
║     return BestMotifs                                                          ║
╚════════════════════════════════════════════════════════════════════════════════╝

Implement GreedyMotifSearch
Given:
Integers k and t, followed by a collection of strings Dna.

Return:
A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t).
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

Plan 1.
  * implementing GreedyMotifSearch using "Profile-most probable k-mer" function in BA2C
  * Pseudocode:
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

  * sample dataset
    3 5
    GGCGTTCAGGCA
    AAGAATCAGTCA
    CAAGGAGTTCGC
    CACGTCAATCAC
    CAATAATATTCG

═════════════════════════════════════════════════

References:
- How do I zip two arrays in JavaScript?
  https://stackoverflow.com/questions/22015684/how-do-i-zip-two-arrays-in-javascript
- Array.prototype.map()
  https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/map#parameters
- Javascript max() for array column
  https://stackoverflow.com/questions/11190407/javascript-max-for-array-column
*/

// BA03A: all k-mers, generator version
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

// Profile(Motifs)
function get_profile(motifs) {
    /**
     * [str] -> {char:[float]}
     * returns a dictionary (object) of profile
     * embedding 'COUNT(Motifs)' function
     * >>> get_profile(["ACGTC", "ATTTT"])
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
    for (const nuc in dic_profile) {
        dic_profile[nuc] = dic_profile[nuc].map(num => num / t);
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

// BA02D: greedy motif search
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
            const profile = get_profile(motifs);
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
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02d.txt').toString().split("\n");
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
