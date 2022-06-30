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
  ╔═════════════════════════════════════════════════════════════════════════════════════════╗
  ║ GIBBSSAMPLER(Dna, k, t, N)                                                              ║
  ║     randomly select k-mers Motifs = (Motif[1], ..., Motif[t]) in each string from Dna   ║
  ║     BestMotifs <- Motifs                                                                ║
  ║     for j <- 1 to N                                                                     ║
  ║         i <- Random(t)                                                                  ║
  ║         Profile <- profile matrix constructed from all strings in Motifs                ║
  ║                    except for Motif[i]                                                  ║
  ║         Motif[i] <- Profile-randomly generated k-mer in the i-th sequence               ║
  ║         if Score(Motifs) < Score(BestMotifs)                                            ║
  ║             BestMotifs <- Motifs                                                        ║
  ║     return BestMotifs                                                                   ║
  ╚═════════════════════════════════════════════════════════════════════════════════════════╝

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
- How to select all other values in an array except the ith element?
  https://stackoverflow.com/questions/15361189/how-to-select-all-other-values-in-an-array-except-the-ith-element
- How to choose a "weighted random" array element in Javascript?
  https://stackoverflow.com/questions/43566019/how-to-choose-a-weighted-random-array-element-in-javascript
  Python  ==> random.choices(kmers, weights=norm_probs, k=1)
*/

// BA03A: all k-mers
function* all_kmers_gen(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length - k + 1; i++) yield text.slice(i, i + k);
}

// BA03A: all k-mers
function all_kmers(text, k) {
    /**
     * (str,int) -> [str]
     * returns a list of kmers
     */
    let kmers = [];
    for (let i = 0; i < text.length - k + 1; i++) kmers.push(text.slice(i, i + k));
    return kmers;
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

// Profile, with pseudocount
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

// BA02G: profile randomly generated k-mer, TODO: fix this
function profile_randomly_generated_kmer(dna, k, profile) {
    /**
     * (str,int,{char:[float]})
     * returns a profile randomly generated k-mer
     */
    // get array of probability weights
    let kmers = [];
    let weights = [];
    var i;
    for (i = 0; i < dna.length - k + 1; i++) {
        const kmer = dna.slice(i, i + k);
        kmers.push(kmer);
        weights.push(get_probability(kmer, profile));
    }
    // cummuative weight slice
    for (i = 0; i < kmers.length; i++) {
        weights[i] += weights[i - 1] || 0;
    }
    // get an element randomly by probability weights
    let rand = Math.random() * weights[weights.length - 1];
    for (i = 0; i < weights.length; i++) {
        if (weights[i] > rand) break;
    }
    return kmers[i];
}

// BA02G: Gibbs sampler
function gibbs_sampler(dnas, k, t, n) {
    /**
     * ([str],int,int,int) -> [str]
     * rturns gibbs sampler
     */
    let motifs = [];
    for (const dna of dnas) {
        const idx = Math.floor(Math.random() * (dna.length - k + 1));
        motifs.push(dna.slice(idx, idx + k));
    }
    let best_motifs = motifs;
    for (const x of Array(n).keys()) {  // range(n) in Python
        const idx = Math.floor(Math.random() * t);
        temp_motifs = motifs.filter((value, index) => index !== idx);
        const profile = get_profile_pseudocount(temp_motifs);
        motifs[idx] = profile_randomly_generated_kmer(dnas[idx], k, profile);
        if (get_score(motifs) < get_score(best_motifs)) best_motifs = motifs;
    }
    return best_motifs;
}

// repeat 20 times
function repeat_n_times(dnas, k, t, n, repeat) {
    let best_motifs = [];
    let best_score = Infinity;
    for (i = 0; i < repeat; i++) {
        let motifs = gibbs_sampler(dnas, k, t, n);
        let score = get_score(motifs);
        if (score < best_score) {
            best_score = score;
            best_motifs = motifs;
        }
    }
    return best_motifs;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02g.txt').toString().split("\n");
    const [k, t, n] = lines[0].split(/\s/).map(Number);
    const dnas = lines.slice(1);

    const startTime = performance.now()
    const motifs = repeat_n_times(dnas, k, t, n, 20);
    for (const motif of motifs) {
        console.log(motif);
    }
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
