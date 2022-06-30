/************************************************
Rosalind: BA1I
Find the Most Frequent Words with Mismatches in a String

We defined a mismatch in "Compute the Hamming Distance Between Two Strings". We
now generalize "Find the Most Frequent Words in a String" to incorporate
mismatches as well.

Given strings Text and Pattern as well as an integer d, we define Countd(Text,
Pattern) as the total number of occurrences of Pattern in Text with at most d
mismatches. For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because
AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA,
AAACA, and AAAGA. Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in Text is simply a string Pattern
maximizing Count_d(Text, Pattern) among all k-mers. Note that Pattern does not
need to actually appear as a substring of Text; for example, AAAAA is the most
frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though AAAAA
does not appear exactly in this string. Keep this in mind while solving the
following problem.

Frequent Words with Mismatches Problem
Find the most frequent k-mers with mismatches in a string.

Given: A string Text as well as integers k and d.

Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
GATG ATGC ATGT

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

  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
  ║    Neighborhood <- the set consisting of the single string Pattern    ║
  ║    for i = 1 to |Pattern|                                             ║
  ║      symbol <- i-th nucleotide of Pattern                             ║
  ║      for each nucleotide x different from symbol                      ║
  ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
  ║        add Neighbor to Neighborhood                                   ║
  ║    return Neighborhood                                                ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════╗
  ║  NEIGHBORS(Pattern, d)                                   ║
  ║    if d = 0                                              ║
  ║      return {Pattern}                                    ║
  ║    if |Pattern| = 1                                      ║
  ║      return {A, C, G, T}                                 ║
  ║    Neighborhood <- an empty set                          ║
  ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
  ║    for each string Text from SuffixNeighbors             ║
  ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
  ║        for each nucleotide x                             ║
  ║          add x + Text to Neighborhood                    ║
  ║      else                                                ║
  ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
  ║    return Neighborhood                                   ║
  ╚══════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════╗
  ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
  ║    Neighborhood <- set consisting of single string Pattern   ║
  ║    for j = 1 to d                                            ║
  ║      for each string Pattern' in Neighborhood                ║
  ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
  ║        remove duplicates from Neighborhood                   ║
  ║    return Neighborhood                                       ║
  ╚══════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  FREQUENTWORDSWITHMISMATCHES(Text, k, d)                              ║
  ║    FrequentPatterns an empty set                                      ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLOSE(i) <- 0                                                    ║
  ║      FREQUENCYARRAY <- 0                                              ║
  ║    for i <- 0 to |Text| - k                                           ║
  ║      Neighborhood <- NEIGHBORS(Text(i, k), d)                         ║
  ║      for each Pattern from Neighborhood                               ║
  ║        index <- p2n(Pattern)                                          ║
  ║        CLOSE(index) <- 1                                              ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLOSE(i) = 1                                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        FREQUENCYARRAY(i) <- APPROXIMATEPATTERNCOUNT(Text, Pattern, d) ║
  ║    maxCount maximal value in FREQUENCYARRAY                           ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if FREQUENCYARRAY(i) = maxCount                                  ║
  ║        Pattern <- n2p(i, k)                                           ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════════════════╗
  ║  FINDINGFREQUENTWORDSWITHMISMATCHESBYSORTING(Text, k, d)                 ║
  ║    FrequentPatterns <- an empty set                                      ║
  ║    Neighborhoods <- an empty list                                        ║
  ║    for i <- 0 to |Text| - k                                              ║
  ║      add NEIGHBORS(Text(i, k), d) to Neighborhoods                       ║
  ║    form an array NEIGHBORHOODARRAY holding all strings in Neighborhoods  ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      Pattern <- NEIGHBORHOODARRAY(i)                                     ║
  ║      INDEX(i) <-  p2n(Pattern)                                           ║
  ║      COUNT(i) <- 1                                                       ║
  ║    SORTEDINDEX SORT(INDEX)                                               ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║      if SORTEDINDEX(i) = SORTEDINDEX(i + 1)                              ║
  ║        COUNT(i + 1) <- COUNT(i) + 1                                      ║
  ║    maxCount maximum value in array COUNT                                 ║
  ║    for i <- 0 to |Neighborhoods| - 1                                     ║
  ║     if COUNT(i) = maxCount                                               ║
  ║       Pattern <- n2p(SORTEDINDEX(i), k)                                  ║
  ║       add Pattern to FrequentPatterns                                    ║
  ║    return FrequentPatterns                                               ║
  ╚══════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════
References:
- Ignoring values in Go's range
  https://stackoverflow.com/questions/49191925/ignoring-values-in-gos-range

************************************************/

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

// helper fx: all approximate occurrences
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

// BA01I: frequent words mismatches
func FreqWordMismatches(text string, k, d int) []string {
	freqWordsMap := make(map[string]bool) // no Set in Golang
	var freqWords []string
	close := make(map[int]int)
	freqMap := make(map[int]int)
	for i := 0; i < len(text)-k+1; i++ {
		neighborhood := Neighbors(text[i:i+k], d)
		for _, pat := range neighborhood {
			num := PatternToNumber(pat)
			close[num] = 1
		}
	}
	maxCount := 0
	for i := range close { // simplify range expression
		pattern := NumberToPattern(i, k)
		freqMap[i] = AllAproximatePatternCount(text, pattern, d)
		if freqMap[i] > maxCount {
			maxCount = freqMap[i]
		}
	}
	for i, val := range freqMap {
		if val == maxCount {
			pattern := NumberToPattern(i, k)
			freqWordsMap[pattern] = true
		}
	}
	for pat := range freqWordsMap { //simplify range expression
		freqWords = append(freqWords, pat)
	}
	return freqWords
}

// code by sonia
func FrequentWordsWithMismatches(Text string, k, d int) (m []string) {
	c := map[string][]string{}
	f := map[string]int{}
	max := 0
	for i, j := 0, k; j <= len(Text); i, j = i+1, j+1 {
		k0 := Text[i:j]
		v, ok := c[k0]
		if !ok {
			v = Neighbors(k0, d)
			c[k0] = v
		}
		for _, kmer := range v {
			n := f[kmer] + 1
			f[kmer] = n
			switch {
			case n > max:
				m = []string{kmer}
				max = n
			case n == max:
				m = append(m, kmer)
			}
		}
	}
	return
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01i.txt")
	text := lines[0]
	var nums []int
	for _, val := range strings.Split(lines[1], " ") {
		num, _ := strconv.Atoi(val)
		nums = append(nums, num)
	}
	k := nums[0]
	d := nums[1]

	// code by sonia - very fast
	start := time.Now()
	for _, kmer := range FrequentWordsWithMismatches(text, k, d) {
		fmt.Printf("%s ", kmer)
	}
	fmt.Println()
	elapsed := time.Since(start)
	fmt.Println("Execution time:", elapsed) // 122.4585 ms

	// slower than sonia's code - Why?
	start = time.Now()
	for _, kmer := range FreqWordMismatches(text, k, d) {
		fmt.Printf("%s ", kmer)
	}
	fmt.Println()
	elapsed = time.Since(start)
	fmt.Println("Execution time:", elapsed) // 725.9799 ms
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

// BA01N: Neighbors
func Neighbors(pattern string, distance int) []string {
	if distance == 0 {
		return []string{pattern}
	}
	if len(pattern) == 1 {
		return []string{"A", "C", "G", "T"}
	}
	neighborhood := []string{}
	suffixPattern := pattern[1:]
	suffixNeighbors := Neighbors(suffixPattern, distance)
	for _, neighbor := range suffixNeighbors {
		if HammingDistance(suffixPattern, neighbor) < distance {
			for _, nuc := range []string{"A", "C", "G", "T"} {
				neighborhood = append(neighborhood, string(nuc)+neighbor)
			}
		} else {
			neighborhood = append(neighborhood, string(pattern[0])+neighbor)
		}
	}
	// sort.Strings(neighborhood)
	return neighborhood
}

// BA01H: all approximate occurrences
func AllAproximateOccurrences(pattern, genome string, d int) []int {
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

// BA01L: PatternToNumber
func PatternToNumber(pattern string) int {
	var num = 0
	var l = len(pattern)
	for i, nuc := range pattern {
		temp := int(math.Pow(4, float64(l-i-1)))
		if nuc == 'A' {
			num += temp * 0
		} else if nuc == 'C' {
			num += temp * 1
		} else if nuc == 'G' {
			num += temp * 2
		} else {
			num += temp * 3
		}
	}
	return num
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
