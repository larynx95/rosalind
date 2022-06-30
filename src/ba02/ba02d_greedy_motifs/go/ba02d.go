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

Info:
    Profile-most Probable k-mer Problem
    Find a Profile-most probable k-mer in a string.
    (profile-most probable k-mer is not consensus.)

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

═════════════════════════════════════════════════
References:
- How to print struct variables in console?
  https://stackoverflow.com/questions/24512112/how-to-print-struct-variables-in-console
  fmt.Printf("%+v\n", yourProject)

*/

package main

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

// implementing 'Pr(pattern, profile)' function
func GetProbability(pattern string, profile map[rune][]float64) float64 {
	pr := 1.0
	for i, nuc := range pattern {
		pr *= profile[nuc][i]
	}
	return pr
}

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
func GetProfile(motifs []string) map[rune][]float64 {
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
			profile[nuc] = append(profile[nuc], float64(temp[nuc])/float64(t))
		}
	}
	return profile
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

// BA02D: greedy motif search
func GreedyMotifSearch(dnas []string, k int, t int) []string {
	// first k-mer in each dna string
	var bestMotifs []string
	for i := 0; i < t; i++ {
		bestMotifs = append(bestMotifs, dnas[i][:k])
	}
	//
	for i := 0; i < len(dnas[0])-k+1; i++ {
		motif_1 := dnas[0][i : i+k]
		motifs := []string{motif_1}
		for j := 1; j < t; j++ {
			profile := GetProfile(motifs)
			motif_i := ProfileMostProbableKmer(dnas[j], k, profile)
			motifs = append(motifs, motif_i)
		}
		if GetScore(motifs) < GetScore(bestMotifs) {
			bestMotifs = motifs
		}
	}
	return bestMotifs
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02d.txt")
	temp := make([]int, 2)
	for i, str := range strings.Split(lines[0], " ") {
		temp[i], _ = strconv.Atoi(str)
	}
	k := temp[0]
	t := temp[1]
	dnas := lines[1:]

	start := time.Now()
	for _, motif := range GreedyMotifSearch(dnas, k, t) {
		println(motif)
	}
	println("Execution time:", time.Since(start))
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
