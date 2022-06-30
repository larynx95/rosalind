/*
Rosalind: BA2C
Find a Profile-most Probable k-mer in a String

Given a profile matrix Profile,
we can evaluate the probability of every k-mer in a string Text
and find a Profile-most probable k-mer in Text,
i.e., a k-mer that was most likely to have been generated
by Profile among all k-mers in Text.
For example, ACGGGGATTACC is the Profile-most probable 12-mer in GGTACGGGGATTACCT.
Indeed, every other 12-mer in this string has probability 0.

In general, if there are multiple Profile-most probable k-mers in Text,
then we select the first such k-mer occurring in Text.

Given: A string Text, an integer k, and a 4 * k matrix Profile.

Return: A Profile-most probable k-mer in Text.
(If multiple answers exist, you may return any one.)

Sample Dataset
ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
5
0.2 0.2 0.3 0.2 0.3
0.4 0.3 0.1 0.5 0.1
0.3 0.3 0.5 0.2 0.4
0.1 0.2 0.1 0.1 0.2

Sample Output
CCGAG

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
      -> PREV: Median string problem (BA2B) - not fast enough... find another algorithm
      ↓
    * GreedyMotifSearch
      -> HERE: Profile-most probable k-mer problem (BA2C)
      -> NEXT: GreedyMotifSearch (BA2D)

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

  * What is the "Greedy algorithm"?
    ; A greedy algorithm is an approach for solving a problem
      by selecting the best option available at the moment.
      It doesn't worry whether the current best result will bring the overall optimal result.
    ; Profile-most probable k-mer finding is a kind of greedy algorithm.

  * return type? Count(Motifs), Profile(Motifs)

Plan 1.
  * steps
    a. get all k-mers from the given text
    b. get probability score for each kmer
    c. find kmers with

═════════════════════════════════════════════════

References:
- 2D slice
  - Slice of Slices in Golang
    https://www.geeksforgeeks.org/slice-of-slices-in-golang/
    a := make([][]int , 3)
  - What is a concise way to create a 2D slice in Go?
    https://stackoverflow.com/questions/39804861/what-is-a-concise-way-to-create-a-2d-slice-in-go
  - How to append to a 2d slice
    https://stackoverflow.com/questions/52706495/how-to-append-to-a-2d-slice
- Parsing float number
  - func ParseFloat
    https://pkg.go.dev/strconv#ParseFloat

*/

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

// implementing 'Pr(pattern, profile)' function
func PrScore(pattern string, profile map[rune][]float64) float64 {
	pr := 1.0
	for i, nuc := range pattern {
		pr *= profile[nuc][i]
	}
	return pr
}

// BA02C: profile most probable k-mer
func ProfileMostProbableKmer(dna string, k int, profile map[rune][]float64) string {
	var answer string
	var maxScore float64 = math.Inf(-1)
	for i := 0; i < len(dna)-k+1; i++ {
		kmer := dna[i : i+k]
		score := PrScore(kmer, profile)
		if score > maxScore {
			maxScore = score
			answer = kmer
		}
	}
	return answer
}

// main
func main() {
	// readlines
	lines := ReadLines("/home/wsl/rosalind/data/ba02c.txt")
	text := lines[0]
	k, _ := strconv.Atoi(lines[1])
	profile := make(map[rune][]float64)
	nucleotides := [4]rune{'A', 'C', 'G', 'T'}
	for i, line := range lines[2:] {
		var values []float64
		for _, str := range strings.Split(line, " ") {
			f, _ := strconv.ParseFloat(str, 64)
			values = append(values, f)
		}
		profile[nucleotides[i]] = values
	}

	// execution time
	start := time.Now()
	println(ProfileMostProbableKmer(text, k, profile))
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

// get all k-mers
func AllKmers(text string, k int) []string {
	var kmers []string
	for i := 0; i < len(text)-k+1; i++ {
		kmers = append(kmers, text[i:i+k])
	}
	return kmers
}
