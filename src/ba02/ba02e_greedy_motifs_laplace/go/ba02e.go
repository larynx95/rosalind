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
func GreedyMotifSearchPseudocount(dnas []string, k int, t int) []string {
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
			profile := GetProfile(motifs, 1)
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
	lines := ReadLines("/home/wsl/rosalind/data/ba02e.txt")
	temp := make([]int, 2)
	for i, str := range strings.Split(lines[0], " ") {
		temp[i], _ = strconv.Atoi(str)
	}
	k := temp[0]
	t := temp[1]
	dnas := lines[1:]

	start := time.Now()
	for _, motif := range GreedyMotifSearchPseudocount(dnas, k, t) {
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
