/************************************************
Rosalind: BA1J
Find Frequent Words with Mismatches and Reverse Complements

We now extend “Find the Most Frequent Words with Mismatches in a String” to find
frequent words with both mismatches and reverse complements. Recall that Pattern
refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
Find the most frequent k-mers (with mismatches and reverse complements) in a DNA
string.

Given:
A DNA string Text as well as integers k and d.

Return:
All k-mers Pattern maximizing the sum
Count_d(Text, Pattern) + Count_d(Text, rPattern) over all possible k-mers.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
ATGT ACAT

═════════════════════════════════════════════════
References:
- named return value
  Empty return in func with return value in golang [duplicate]
  https://stackoverflow.com/questions/45239409/empty-return-in-func-with-return-value-in-golang

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

/* code by sonia */

func RevC(Pattern string) string {
	rc := make([]byte, len(Pattern))
	rcx := len(rc)
	for _, b := range []byte(Pattern) {
		rcx--
		rc[rcx] = ^b&2>>1*17 | 4 ^ b
	}
	return string(rc)
}

func Neighbors1(s string, h int) []string {
	v := []string{s}        // result list
	const sym = "A C T G"   // used for encoding symbols from indexes
	buf := []byte(s)        // buffer for generating neighbors
	var f func([]byte, int) // recursive function
	f = func(t []byte, h int) {
		for i := 0; i < len(t); i++ {
			sub := t[i:]
			b := sub[0]
			for vb := byte(0); ; vb += 2 {
				if vb == b&6 {
					vb += 2
				}
				if vb == 8 {
					break
				}
				sub[0] = sym[vb]
				v = append(v, string(buf))
				if h > 1 && len(sub) > 1 {
					f(sub[1:], h-1)
				}
			}
			sub[0] = b
		}
	}
	f(buf, h)
	return v
}

func FrequentWordsWithMismatchesRC(Text string, k, d int) (m []string) {
	c := map[string][]string{}
	f := map[string]int{}
	max := 0
	tally := func(k0 string) {
		v, ok := c[k0]
		if !ok {
			v = Neighbors1(k0, d)
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
	revc := RevC(Text)
	for i, j := 0, k; j <= len(Text); i, j = i+1, j+1 {
		tally(Text[i:j])
		tally(revc[i:j])
	}
	return
}

/* my code */

// BA01J: FrequentWordsMismatchesReverseComplement
func FreqWordMismatchesRC(text string, k, d int) []string {
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
		freqMap[i] = AllAproximatePatternCount(text, pattern, d) + AllAproximatePatternCount(text, ReverseComplement(pattern), d)
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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01j.txt")
	text := lines[0]
	var nums []int
	for _, val := range strings.Split(lines[1], " ") {
		num, _ := strconv.Atoi(val)
		nums = append(nums, num)
	}
	k := nums[0]
	d := nums[1]

	// check time
	start := time.Now()
	for _, kmer := range FrequentWordsWithMismatchesRC(text, k, d) {
		fmt.Printf("%s ", kmer)
	}
	fmt.Println()
	elapsed := time.Since(start)
	fmt.Println("Execution time:", elapsed)

	// check time - slower
	start = time.Now()
	for _, kmer := range FreqWordMismatchesRC(text, k, d) {
		fmt.Printf("%s ", kmer)
	}
	fmt.Println()
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

// helper fx
func ReverseString(pattern string) string {
	var reversed string
	for i := len(pattern) - 1; i >= 0; i-- {
		reversed = reversed + string(pattern[i])
	}
	return reversed
}

// BA01C: ReverseComplement
func ReverseComplement(pattern string) string {
	var complement_string string
	reversed := ReverseString(pattern)
	for i := 0; i < len(reversed); i++ {
		var nucleotide byte
		switch reversed[i] {
		case 'A':
			nucleotide = 'T'
		case 'C':
			nucleotide = 'G'
		case 'G':
			nucleotide = 'C'
		case 'T':
			nucleotide = 'A'
		}
		complement_string += string(nucleotide)
	}
	return complement_string
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

// helper fx
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
