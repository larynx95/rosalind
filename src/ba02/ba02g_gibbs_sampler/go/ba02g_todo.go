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
- rand package
  https://go.dev/play/p/RaoKyUd9tgC
- Weighted Random: algorithms for sampling from discrete probability distributions
  https://zliu.org/post/weighted-random/
- Fast weighted random selection for Go
  https://golangexample.com/fast-weighted-random-selection-for-go/

*/

package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

const PSEUDOCOUNT = 1
const MaxUint = ^uint(0)
const MinUint = 0
const MaxInt = int(MaxUint >> 1)
const MinInt = -MaxInt - 1

// Count fx
func GetCount(motifs []string) map[rune][]int {
	cnt := map[rune][]int{'A': {}, 'C': {}, 'G': {}, 'T': {}}
	for i := 0; i < len(motifs[0]); i++ {
		dic := map[rune]int{'A': 0, 'C': 0, 'G': 0, 'T': 0}
		for j := 0; j < len(motifs); j++ {
			if motifs[j][i] == 'A' {
				dic['A'] += 1
			} else if motifs[j][i] == 'C' {
				dic['C'] += 1
			} else if motifs[j][i] == 'G' {
				dic['G'] += 1
			} else if motifs[j][i] == 'T' {
				dic['T'] += 1
			}
		}
		for _, nuc := range "ACGT" {
			cnt[nuc] = append(cnt[nuc], dic[nuc])
		}
	}
	return cnt
}

// Score
func GetScore(motifs []string) int {
	t := len(motifs)
	n := len(motifs[0])
	countMap := GetCount(motifs)
	score := 0
	for i := 0; i < n; i++ {
		maxCount := 0
		for _, nuc := range "ACGT" {
			cnt := countMap[nuc][i]
			if cnt > maxCount {
				maxCount = cnt
			}
		}
		score += t - maxCount
	}
	return score
}

// Profile
func GetProfile(motifs []string, pseudocount int) map[rune][]float64 {
	t := len(motifs)
	l := len(motifs[0])
	profile := map[rune][]float64{'A': {}, 'C': {}, 'G': {}, 'T': {}}
	for i := 0; i < l; i++ {
		temp := map[rune]int{'A': 0, 'C': 0, 'G': 0, 'T': 0}
		for j := 0; j < t; j++ {
			if motifs[j][i] == 'A' {
				temp['A'] += 1
			} else if motifs[j][i] == 'C' {
				temp['C'] += 1
			} else if motifs[j][i] == 'G' {
				temp['G'] += 1
			} else if motifs[j][i] == 'T' {
				temp['T'] += 1
			}
		}
		for _, nuc := range "ACGT" {
			profile[nuc] = append(profile[nuc], float64(temp[nuc]+pseudocount)/float64(t+t*pseudocount))
		}
	}
	return profile
}

// implementing 'Pr(pattern, profile)' function
func GetProbability(pattern string, profile map[rune][]float64) float64 {
	pr := 1.0
	for i, nuc := range pattern {
		pr *= profile[nuc][i]
	}
	return pr
}

// BA02G: profile randomly generated k-mer
func ProfileRandomlyGeneratedKmer(dna string, k int, profile map[rune][]float64) string {
	// create two slices for k-mers and their weights
	kmers := []string{}
	weights := []float64{}
	for i := 0; i < len(dna)-k+1; i++ {
		kmer := dna[i : i+k]
		kmers = append(kmers, kmer)
		weights = append(weights, GetProbability(kmer, profile))
	}
	// get cumulative weight slice
	for i := 1; i < len(weights); i++ {
		weights[i] += weights[i-1]
	}
	// TODO: What is goint on?
	rnum := rand.Float64() * weights[len(weights)-1]
	for i := 0; i < len(weights); i++ {
		if weights[i] > rnum {
			return kmers[i]
		}
	}
	return kmers[len(kmers)-1]
}

// BA02G: Gibbs Sampler, TODO: Fix this function
func GibbsSampler(dnas []string, k int, t int, n int) []string {
	// randomly selected k-mers Motifs in each dna string
	motifs := []string{}
	for _, dna := range dnas {
		rIdx := rand.Intn(len(dna) - k + 1)
		motifs = append(motifs, dna[rIdx:rIdx+k])
	}
	bestMotifs := motifs
	// repeat n times
	for j := 0; j < n; j++ {
		// get random integer rIdx from range [0..t)
		rIdx := rand.Intn(t)
		// remove Motifs[i] from Motifs
		motifsRemoved := []string{}
		for i, motif := range motifs {
			if i != rIdx {
				motifsRemoved = append(motifsRemoved, motif)
			}
		}
		// get profile from all Motifs except for Motifs[rIdx]
		profile := GetProfile(motifsRemoved, PSEUDOCOUNT)
		motifs[rIdx] = ProfileRandomlyGeneratedKmer(dnas[rIdx], k, profile)
		if GetScore(motifs) < GetScore(bestMotifs) {
			bestMotifs = motifs
		} else {
			return bestMotifs
		}
	}
	return bestMotifs
}

// repeat n times
func repeat_n_times(dnas []string, k int, t int, n int, repeat int) []string {
	best_score := MaxInt
	best_motifs := []string{}
	for i := 0; i < repeat; i++ {
		motifs := GibbsSampler(dnas, k, t, n)
		score := GetScore(motifs)
		if score < best_score {
			best_score = score
			best_motifs = motifs
		}
	}
	return best_motifs
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02g.txt")
	temp := make([]int, 3)
	for i, str := range strings.Split(lines[0], " ") {
		temp[i], _ = strconv.Atoi(str)
	}
	k := temp[0]
	t := temp[1]
	n := temp[2]
	dnas := lines[1:]

	// check time
	start := time.Now()
	for _, motif := range repeat_n_times(dnas, k, t, n, 20) {
		println(motif)
	}
	fmt.Println("Execution time:", time.Since(start))

}

// helper fx: read lines
func ReadLines(path string) []string {
	// open file, close file
	file, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer file.Close()

	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	return lines
}

// ═════════════════════════════════════════════════
// weighted random
// https://zliu.org/post/weighted-random/

func WeightedRandomS1(weights []float32) int {
	if len(weights) == 0 {
		return 0
	}
	var choices []int
	for i, w := range weights {
		wi := int(w * 10)
		for j := 0; j < wi; j++ {
			choices = append(choices, i)
		}
	}
	return choices[rand.Int()%len(choices)]
}

func WeightedRandomS2(weights []float32) int {
	if len(weights) == 0 {
		return 0
	}
	var sum float32 = 0.0
	for _, w := range weights {
		sum += w
	}
	r := rand.Float32() * sum
	for i, w := range weights {
		r -= w
		if r < 0 {
			return i
		}
	}
	return len(weights) - 1
}

func WeightedRandomS3(weights []float32) int {
	n := len(weights)
	if n == 0 {
		return 0
	}
	cdf := make([]float32, n)
	var sum float32 = 0.0
	for i, w := range weights {
		if i > 0 {
			cdf[i] = cdf[i-1] + w
		} else {
			cdf[i] = w
		}
		sum += w
	}
	r := rand.Float32() * sum
	var l, h int = 0, n - 1
	for l <= h {
		m := l + (h-l)/2
		if r <= cdf[m] {
			if m == 0 || (m > 0 && r > cdf[m-1]) {
				return m
			}
			h = m - 1
		} else {
			l = m + 1
		}
	}
	return -1
}
