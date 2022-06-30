/*
Rosalind: BA2B
Find a Median String

Given a k-mer Pattern and a longer string Text,
we use d(Pattern, Text) to denote the minimum Hamming distance
between Pattern and any k-mer in Text,

d(Pattern, Text) = min                 HammingDistance(Pattern, Pattern')
                  all k-mers Pattern' in Text

    example:
    d(GATTCTCA, gcaaaGACGCTGAccaa) = 3

Given a k-mer Pattern and a set of strings
Dna = {Dna_1, ... , Dna_t},
we define d(Pattern, Dna) as the sum of distances
between Pattern and all strings in Dna,

                    t
    d(Pattern, Dna) = Σ  d(Pattern, Dna_i)
                    i=0

    example:
           d(AAA, Dna) = 1 + 1 + 2 + 0 + 1 = 5
           ttaccttAAC   1
           gATAtctgtc   1
    Dna    ACGgcgttcg   2
           ccctAAAgag   0
           cgtcAGAggt   1

Our goal is to find a k-mer Pattern
that minimizes d(Pattern, Dna) over all k-mers Pattern,
the same task that the Equivalent Motif Finding Problem is trying to achieve.
We call such a k-mer a "median string" for Dna.

Median String Problem
Find a median string.

Given: An integer k and a collection of strings Dna.

Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern.
(If multiple answers exist, you may return any one.)

Sample Dataset
3
AAATTGACGCAT
GACGACCACGTT
CGTCAGCGCCTG
GCTGAGCACCGG
AGTACGGGACAG

Sample Output
GAC

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> PREV: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference? (very important)
      -> HERE: "Median string problem" (BA2B) - not fast enough... find better algorithm
      ↓
    * GreedyMotifSearch
      -> NEXT: Profile-most probable k-mer problem (BA2C)

Info:
- Counting row-by-row returns exactly the same result as counting col-by-col.
                                                               (2)
                          T  C  G  G  G  G  g  T  T  T  t  t   3 row-by-row
                          c  C  G  G  t  G  A  c  T  T  a  C   4
                          a  C  G  G  G  G  A  T  T  T  t  C   2
                          T  t  G  G  G  G  A  c  T  T  t  t   4
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C   3
                          T  t  G  G  G  G  A  c  T  T  C  C   2
                          T  C  G  G  G  G  A  T  T  c  a  t   3
                          T  C  G  G  G  G  A  T  T  c  C  t   2
                          T  a  G  G  G  G  A  a  c  T  a  C   4
                          T  C  G  G  G  t  A  T  a  a  C  C + 3   Interesting!
    SCORE(Motifs)   (1)   3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30  result score is the same!
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C
                          col-by-col
                  ╔══════════════════════════════════════════════════╗
                  ║   SCORE(Motifs) == d(CONSENSUS(Motifs), Motifs)  ║
                  ║   col-by-col       row-by-row                    ║
                  ╚══════════════════════════════════════════════════╝
   (1) given motifs -> consensus
      ; Motifs -> SCORE(Motifs) -> CONSENSUS(Motifs)
      ; Motif Finding Problem
   (2) condidates of consensus -> best possible collection Motifs for each consensus string
      ; CONSENSUS(Motifs) -> d(CONSENSUS(Motifs), Motifs) -> Motifs
      ; Equivalent Motif Finding Problem

                      (1) first approach             (2) second approach
                      T C G G G G g T T T t t        T C G G G G g T T T t t   3
                      c C G G t G A c T T a C        c C G G t G A c T T a C   4
                      a C G G G G A T T T t C        a C G G G G A T T T t C   2
                      T t G G G G A c T T t t        T t G G G G A c T T t t   4
    Motifs            a a G G G G A c T T C C        a a G G G G A c T T C C   3
                      T t G G G G A c T T C C        T t G G G G A c T T C C   2
                      T C G G G G A T T c a t        T C G G G G A T T c a t   3
                      T C G G G G A T T c C t        T C G G G G A T T c C t   2
                      T a G G G G A a c T a C        T a G G G G A a c T a C   4
                      T C G G G t A T a a C C        T C G G G t A T a a C C + 3
    SCORE(Motifs)     3+4+0+0+1+1+1+5+2+3+6+4 = 30   3+4+0+0+1+1+1+5+2+3+6+4= 30

- comparing three problems
  ┌──────────────────────────────────────────────────────────────────┐
  │ Motif Finding Problem:                                           │ O(n^t * k * t)
  │   Given a collection of strings,                                 │
  │   find a set of k-mers, one from each string,                    │
  │   that minimizes the score of the resulting motif.               │
  │                                                                  │
  │ Input:                                                           │
  │   A collection of strings Dna and an integer k.                  │
  │ Output:                                                          │
  │   A collection "Motifs" of k-mers, one from each string in Dna,  │ Motifs
  │   minimizing SCORE(Motifs) among all possible choices of k-mers. │
  └──────────────────────────────────────────────────────────────────┘
  ┌────────────────────────────────────────────────────────────────────┐
  │ Equivalent Motif Finding Problem:                                  │
  │   Given a collection of strings,                                   │
  │   find a pattern and a collection of k-mers (one from each string) │
  │   that minimizes the distance between all possible patterns        │
  │   and all possible collections of k-mers.                          │
  │                                                                    │
  │ Input:                                                             │
  │   A collection of strings Dna and an integer k.                    │
  │ Output:                                                            │
  │   A k-mer "Pattern" and a collection of k-mers "Motifs",           │ Pattern, Motifs
  │   one from each string in Dna, minimizing d(Pattern, Motifs)       │
  │   among all possible choices of Pattern and Motifs.                │
  └────────────────────────────────────────────────────────────────────┘
  ┌──────────────────────────────────────────────────┐
  │ Median String Problem:                           │
  │   Find a median string.                          │
  │ Input:                                           │
  │   A collection of strings Dna and an integer k.  │
  │ Output:                                          │
  │   A k-mer Pattern minimizing d(Pattern, Dna)     │
  │   among all k-mers Pattern.                      │
  └──────────────────────────────────────────────────┘

- What is the meaning of above comparison?
  ; The col-by-col approach is exactly the same as "Motif finding problem", nothing else.
  ; But row-by-row approach is something different.
    The key observation for solving this Equivalent Motif Finding Problem is that,
    given Pattern, we don't need to explore all possible collections (4^k) Motifs
    in order to minimize d(Pattern, Motifs).

    MOTIFS(AAA, Dna):
                         ttaccttAAC
                         gATAtctgtc
                     Dna ACGgcgttcg
                         ccctAAAgag
                         cgtcAGAggt

Plan 1.
- What is the "median string"?
  ; a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern

- What are the meanings of "d(pattern, text)" and "d(pattern, Dna)"?
  d(pattern, one Dna string)       = x   (minimum hamming distance value)
  d(pattern, multiple Dna strings) = x + y + ... + z

Plan 2.
- pseudocode:
  ╔═══════════════════════════════════════════════════════╗
  ║  MEDIANSTRING(Dna, k)                                 ║
  ║      distance <- infinite                             ║
  ║      for each k-mer Pattern from AA...AA to TT...TT   ║
  ║          if distance > d(Pattern, Dna)                ║
  ║              distance <- d(Pattern, Dna)              ║
  ║              Median <- Pattern                        ║
  ║      return Median                                    ║
  ╚═══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
-
*/

// BA03A: all k-mers, generative version
function* all_kmers(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length - k + 1; i++) yield text.slice(i, i + k);
}

// BA01G: hamming distance
function hamming_distance(astr, bstr) {
    /**
     * (str,str) -> int
     * returns a Hamming distance of two strings (BA1G)
     */
    let distance = 0;
    for (let i = 0; i < astr.length; i++)
        if (astr[i] != bstr[i]) distance++;
    return distance;
}

// d(pattern, text)
function d_pattern_text(pattern, text) {
    /**
     * (str,str) -> int
     * returns d(Pattern,Text),
     * the minimum Hamming distance between Pattern and any k-mer in Text
     */
    let min_hdistance = Infinity;
    const k = pattern.length;
    for (const kmer of all_kmers(text, k)) {
        const dist = hamming_distance(pattern, kmer);
        if (dist < min_hdistance) min_hdistance = dist;
    }
    return min_hdistance;
}

// d(pattern, texts)
function d_pattern_texts(pattern, texts) {
    /**
     * (str,[str]) -> int
     * returns d(Pattern,Texts)
     * the sum of distances between Pattern and all strings in Dna
     * >>> d_pattern_texts("AAA", ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT'])
     *     5
     */
    let result = 0;
    for (const text of texts) {
        result += d_pattern_text(pattern, text);
    }
    return result;
}

// BA01M: number to pattern
function number_to_pattern(num, k) {
    /**
     * (int,int) -> str
     * returns a pattern from a number
     * number_to_pattern(15, 2) == "TT"
     */
    pattern = ""
    for (let i = 0; i < k; i++) {
        temp = Math.floor(num / Math.pow(4, (k - i - 1)));
        if (temp == 0) pattern += 'A';
        else if (temp == 1) pattern += 'C';
        else if (temp == 2) pattern += 'G';
        else if (temp == 3) pattern += 'T';
        num -= temp * Math.pow(4, k - i - 1);
    }
    return pattern
}

// BA02B (BA02H): median string
function median_string(texts, k) {
    /**
     * ([str],int) -> str
     * returns a median string
     * >>> median_string_textbook(['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTACGGGACAG'],3)
     *     {'ACG'}
     */
    let mstr = new Set();
    let min_distance = Infinity;
    for (let i = 0; i < Math.pow(4, k); i++) {
        const pattern = number_to_pattern(i, k);
        const dist = d_pattern_texts(pattern, texts);
        if (min_distance > dist) {
            min_distance = dist;
            mstr = new Set([pattern]);
        } else if (min_distance == dist) {
            mstr.add(pattern);
        }
    }
    return mstr;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02b.txt').toString().split("\n");
    const k = parseInt(lines[0]);
    const dnas = lines.slice(1);

    const startTime = performance.now();
    for (const elem of median_string(dnas, k)) {
        console.log(elem);
    }
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
