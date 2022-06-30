/*
Rosalind: BA1H
Find All Approximate Occurrences of a Pattern in a String

We say that a k-mer Pattern appears as a substring of Text with at most d
mismatches if there is some k-mer substring Pattern' of Text having d or fewer
mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our
observation that a DnaA box may appear with slight variations leads to the
following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem
Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return: All starting positions where Pattern appears as a substring of Text with
at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3

Sample Output
6 7 26 27 78

═════════════════════════════════════════════════

Pseudocode:

  ╔═══════════════════════════════════════════════════╗
  ║  APPROXIMATEPATTERNCOUNT(Text, Pattern, d)        ║
  ║    count 0                                        ║
  ║    for i 0 to |Text| - |Pattern|                  ║
  ║      Pattern <- Text(i, |Pattern|)                ║
  ║      if HAMMINGDISTANCE(Pattern, Pattern’) >= d   ║
  ║        count count + 1                            ║
  ║    return count                                   ║
  ╚═══════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
-
*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"time"
)

// BA01H: ApproximatePatternCount
func AllAproximateOccurrences(genome, pattern string, d int) []int {
	var indices []int
	k := len(pattern)
	for i := 0; i < len(genome)-k+1; i++ {
		subs := genome[i : i+k]
		if HammingDistance(subs, pattern) <= d {
			indices = append(indices, i)
		}
	}
	return indices
}

// BA01H: ApproximatePatternCount
func AllAproximatePatternCount(genome, pattern string, d int) int {
	count := 0
	k := len(pattern)
	for i := 0; i < len(genome)-k+1; i++ {
		subs := genome[i : i+k]
		if HammingDistance(subs, pattern) <= d {
			count += 1
		}
	}
	return count
}

// main
func main() {
	// read file
	lines := ReadLines("/home/wsl/rosalind/data/ba01h.txt")
	pattern := lines[0]
	genome := lines[1]
	d, err := strconv.Atoi(lines[2])
	if err != nil {
		fmt.Println(err)
	}

	// check time
	start := time.Now()
	for _, val := range AllAproximateOccurrences(genome, pattern, d) {
		fmt.Print(val, " ")
	}
	fmt.Println()
	elapsed := time.Since(start)
	fmt.Println("Execution time:", elapsed)

	// check time
	start = time.Now()
	fmt.Println(AllAproximatePatternCount(genome, pattern, d))
	elapsed = time.Since(start)
	fmt.Println("Execution time:", elapsed)
}

// helper fx: readlines
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

// BA01G: hammind distance
func HammingDistance(pattern1, pattern2 string) int {
	var dist = 0
	for i, _ := range pattern1 {
		if pattern1[i] != pattern2[i] {
			dist += 1
		}
	}
	return dist
}
