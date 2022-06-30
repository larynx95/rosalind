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
      -> HERE: "Median string problem" (BA2B, BA2H) - not fast enough... find better algorithm
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

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"
)

// Constant definitions
// https://go.dev/play/p/Cc7lCoutegK
const MaxUint = ^uint(0)
const MinUint = 0
const MaxInt = int(^uint(0) >> 1)
const MinInt = -MaxInt - 1

// BA03A: get all k-mers from a text
func AllKmers(Text string, k int) []string {
	c := make([]string, len(Text)-k+1)
	for i, j := 0, k; j <= len(Text); i, j = i+1, j+1 {
		c[i] = Text[i:j]
	}
	return c
}

// BA01G: HammingDistance
func HammingDistance(pattern1, pattern2 string) int {
	var dist = 0
	for i := range pattern1 {
		if pattern1[i] != pattern2[i] {
			dist += 1
		}
	}
	return dist
}

// implementation of 'd(pattern, text)' fx
func DistancePatternText(pattern, text string) int {
	k := len(pattern)
	distance := k
	for i := 0; i < len(text)-k+1; i++ {
		dist := HammingDistance(pattern, text[i:i+k])
		if dist < distance {
			distance = dist
		}
	}
	return distance
}

// implementation of 'd(pattern, dnas)' fx
func DistancePatternTexts(pattern string, dnas []string) int {
	sum := 0
	for _, dna := range dnas {
		sum += DistancePatternText(pattern, dna)

	}
	return sum
}

// BA01M: NumberToPattern
func NumberToPattern(number int, k int) string {
	var pattern []byte
	for i := 0; i < k; i++ {
		temp := number / int(math.Pow(4, float64(k-i-1)))
		if temp == 0 {
			pattern = append(pattern, 'A')
		} else if temp == 1 {
			pattern = append(pattern, 'C')
		} else if temp == 2 {
			pattern = append(pattern, 'G')
		} else if temp == 3 {
			pattern = append(pattern, 'T')
		}
		number -= temp * int(math.Pow(4, float64(k-i-1)))
	}
	return string(pattern)
}

// BA02B: MedianString
func MedianStrings(dnas []string, k int) []string {
	minDistance := 1000                                 // +infinity
	median := make(map[string]int)                      // use Map, cuz no Set in Golang
	for i := 0; i < int(math.Pow(4, float64(k))); i++ { // for all k-mers (inefficient)
		pattern := NumberToPattern(i, k)                   // number to pattern
		newDistance := DistancePatternTexts(pattern, dnas) // get distance
		if minDistance > newDistance {                     // compare distances
			minDistance = newDistance     // update min distance
			median = make(map[string]int) // if min distance changes, init Set(Map)
		}
		median[pattern] = newDistance // only one median string, TODO: improve this part
	}
	result := []string{}
	for key := range median {
		if median[key] == minDistance {
			result = append(result, key)
		}
	}
	return result
}

// BA02H: distance between pattern and strings
func DistanceBwPatternAndStrings(pattern string, dnas []string) int {
	k := len(pattern)
	distance := 0
	for _, dna := range dnas {
		min_distance := MaxInt
		for i := 0; i < len(dna)-k+1; i++ {
			kmer := dna[i : i+k]
			hdistance := HammingDistance(kmer, pattern)
			if min_distance > hdistance {
				min_distance = hdistance
			}
		}
		distance += min_distance
	}
	return distance
}

// BA02H: median string
func MedianStringBA02H(dnas []string, k int) string {
	minDistance := MaxInt
	var median string
	for i := 0; i < int(math.Pow(4, float64(k))); i++ {
		pattern := NumberToPattern(i, k)
		dist := DistanceBwPatternAndStrings(pattern, dnas)
		if minDistance > dist {
			minDistance = dist
			median = pattern
		}
	}
	return median
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02b.txt")
	k, _ := strconv.Atoi(lines[0])
	dnas := lines[1:]

	start := time.Now()
	medians := MedianStrings(dnas, k)
	fmt.Println(medians[0])
	fmt.Println("Execution time:", time.Since(start))

	start = time.Now()
	println(MedianStringBA02H(dnas, k))
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
