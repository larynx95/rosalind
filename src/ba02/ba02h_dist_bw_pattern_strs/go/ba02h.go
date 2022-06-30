/*
Rosalind: BA2H
Implement DistanceBetweenPatternAndStrings

The first potential issue with implementing MedianString from "Find a Median String" is
writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai),
the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}.
This task is achieved by the following pseudocode.

  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

Sample Dataset
AAA
TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

Sample Output
5

═════════════════════════════════════════════════
Pseudocode:
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ MEDIANSTRING(Dna, k)                                                   ║
  ║     distance 1                                                         ║
  ║     for i 0 to 4^k - 1                                                 ║
  ║         pattern <- NumberToPattern(i, k)                               ║
  ║         if distance > DistanceBetweenPaternAndStrings(Pattern, Dna)    ║
  ║             distance <- DistanceBetweenPatternAndStrings(Patern, Dna)  ║
  ║             Median <- Pattern                                          ║
  ║     return Median                                                      ║
  ╚════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════
References:
- maximum value of Integer
  https://go.dev/play/p/Cc7lCoutegK
*/

package main

import (
	"bufio"
	"math"
	"os"
	"strings"
	"time"
)

// Constant definitions
// https://go.dev/play/p/Cc7lCoutegK
const MaxUint = ^uint(0)
const MinUint = 0
const MaxInt = int(^uint(0) >> 1)
const MinInt = -MaxInt - 1

// BA01G: hamming distance
func HammingDistance(pattern1, pattern2 string) int {
	var dist = 0
	for i, _ := range pattern1 {
		if pattern1[i] != pattern2[i] {
			dist += 1
		}
	}
	return dist
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

// BA02H: distance between pattern and strings
func distance_bw_pattern_strings(pattern string, dnas []string) int {
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
func median_string(dnas []string, k int) string {
	minDistance := MaxInt
	var median string
	for i := 0; i < int(math.Pow(4, float64(k))); i++ {
		pattern := NumberToPattern(i, k)
		dist := distance_bw_pattern_strings(pattern, dnas)
		if minDistance > dist {
			minDistance = dist
			median = pattern
		}
	}
	return median
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02h.txt")
	pattern := lines[0]
	dnas := strings.Split(lines[1], " ")
	// println(pattern)
	// fmt.Printf("%+v\n", dnas)

	// check tiem
	start := time.Now()
	println(distance_bw_pattern_strings(pattern, dnas))
	println("Execution time:", time.Since(start))

	// check tiem
	start = time.Now()
	println(median_string(dnas, 3))
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
