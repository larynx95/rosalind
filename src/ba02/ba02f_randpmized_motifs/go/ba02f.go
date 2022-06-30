/*
Rosalind: BA2F
Implement RandomizedMotifSearch

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
  ║          Profile <- Profile(Motifs)         # get Profile                            ║ Profile
  ║          Motifs <- Motifs(Profile, Dna)     # find best Motifs                       ║ Motifs
  ║          if Score(Motifs) < Score(BestMotifs)                                        ║ Score
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
- Why this repeats the same random number?
  https://stackoverflow.com/questions/45753397/why-this-repeats-the-same-random-number
  golang rand.Int(). Why every time same values? [duplicate]
  https://stackoverflow.com/questions/68203678/golang-rand-int-why-every-time-same-values
    rand.Seed(time.Now().UnixNano())
    println(rand.Intn(10))
- Is 'append' function efficient?
  Remove and adding elements to array in GO lang
  https://stackoverflow.com/questions/33834742/remove-and-adding-elements-to-array-in-go-lang
- How do you time a function in Go and return its runtime in milliseconds?
  https://stackoverflow.com/questions/8350609/how-do-you-time-a-function-in-go-and-return-its-runtime-in-milliseconds
*/

package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

const PSEUDOCOUNT = 1

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

// BA02C: profile most probable k-mer
func ProfileMostProbableKmer(text string, k int, profile map[rune][]float64) string {
	answer := ""
	maxProb := math.Inf(-1) // TODO: Don't use 0.0! This line takes me hours.
	for i := 0; i < len(text)-k+1; i++ {
		kmer := text[i : i+k]
		score := GetProbability(kmer, profile)
		if score > maxProb {
			maxProb = score
			answer = kmer
		}
	}
	return answer
}

// implementing 'Motifs(profile, Pna)' function
func GetMotifs(profile map[rune][]float64, k int, dnas []string) []string {
	motifs := []string{}
	for _, dna := range dnas {
		motifs = append(motifs, ProfileMostProbableKmer(dna, k, profile))
	}
	return motifs
}

// BA02F: randomized motif search
func RandomizedMotifSearch(dnas []string, k int, t int) []string {
	motifs := []string{}
	for _, dna := range dnas {
		rand.Seed(time.Now().UnixNano())
		rIdx := rand.Intn(len(dna) - k + 1)
		motifs = append(motifs, dna[rIdx:rIdx+k])
	}
	bestMotifs := motifs
	for {
		profile := GetProfile(motifs, PSEUDOCOUNT)
		motifs = GetMotifs(profile, k, dnas)
		if GetScore(motifs) < GetScore(bestMotifs) {
			bestMotifs = motifs
		} else {
			return bestMotifs
		}
	}
}

// BA02F: repeat n times
func repeat_n_times(dnas []string, k int, t int, n int) []string {
	minScore := math.Inf(1)
	bestMotifs := []string{}
	for i := 0; i < n; i++ {
		motifs := RandomizedMotifSearch(dnas, k, t)
		score := GetScore(motifs)
		if float64(score) < minScore {
			minScore = float64(score)
			bestMotifs = motifs
		}
	}
	return bestMotifs
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02f.txt")
	temp := make([]int, 2)
	for i, str := range strings.Split(lines[0], " ") {
		temp[i], _ = strconv.Atoi(str)
	}
	k := temp[0]
	t := temp[1]
	dnas := lines[1:]

	// practice random integer
	rand.Seed(time.Now().UnixNano())
	println(rand.Intn(10))

	// check time
	start := time.Now()
	for _, motif := range repeat_n_times(dnas, k, t, 1000) {
		println(motif)
	}
	fmt.Println("Execution time:", time.Since(start)) // 2.7881221s
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
