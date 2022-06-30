/*
Rosalind: BA1B
Find the Most Frequent Words in a String

We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text,
Pattern) among all k-mers. For example, "ACTAT" is a most frequent 5-mer in
"ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most frequent 3-mer of
"CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4

Sample Output
CATG GCAT

═════════════════════════════════════════════════

Pseudocode:

╔═════════════════════════════════════════════════╗
║ COMPUTINGFREQUENCIES(Text, k)                   ║
║   for i <- 0 to 4^k- 1                          ║
║     FREQUENCYARRAY(i) <- 0                      ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- Text(i, k)                       ║
║     j <- PATTERNTONUMBER(Pattern)               ║
║     FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1  ║
║   return FREQUENCYARRAY                         ║
╚═════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════╗
║ FREQUENTWORDS(Text, k)                          ║
║   FrequentPatterns an empty set                 ║
║     for i <- 0 to |Text| - k                    ║
║     Pattern <- the k-mer Text(i, k)             ║
║     COUNT(i) <- PATTERNCOUNT(Text, Pattern)     ║
║   maxCount maximum value in array COUNT         ║
║   for i <- 0 to |Text| - k                      ║
║     if COUNT(i) = maxCount                      ║
║       add Text(i, k) to FrequentPatterns        ║
║   remove duplicates from FrequentPatterns       ║
║   return FrequentPatterns                       ║
╚═════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════╗
║ FASTERFREQUENTWORDS(Text , k)                     ║
║   FrequentPatterns <- an empty set                ║
║   FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k) ║
║   maxCount <- maximal value in FREQUENCYARRAY     ║
║   for i <- 0 to 4k - 1                            ║
║     if FREQUENCYARRAY(i) = maxCount               ║
║       Pattern <- NUMBERTOPATTERN(i, k)            ║
║       add Pattern to the set FrequentPatterns     ║
║   return FrequentPatterns                         ║
╚═══════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════╗
║ FINDINGFREQUENTWORDSBYSORTING(Text , k)              ║
║   FrequentPatterns <- an empty set                   ║
║   for i <- 0 to |Text| - k                           ║
║     Pattern <- Text(i, k)                            ║
║     INDEX(i) <- PATTERNTONUMBER(Pattern)             ║
║     COUNT(i) <- 1                                    ║
║   SORTEDINDEX <- SORT(INDEX)                         ║
║   for i <- 1 to |Text| - k                           ║
║     if SORTEDINDEX(i) = SORTEDINDEX(i - 1)           ║
║       COUNT(i) = COUNT(i - 1) + 1                    ║
║   maxCount <- maximum value in the array COUNT       ║
║   for i <- 0 to |Text| - k                           ║
║     if COUNT(i) = maxCount                           ║
║       Pattern <- NUMBERTOPATTERN(SORTEDINDEX(i), k)  ║
║       add Pattern to the set FrequentPatterns        ║
║   return FrequentPatterns                            ║
╚══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

Plan 1.
- create a dictionary (python) to store key-value pair
  but golang doesn't have dictionary
  a. use 'map'
  b. use binary tree with 'structure'

Topics:
- read file line by line, put data into array or vector
- for-range in golang
- error handling
- two more return values in some golang functions
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"time"
)

// BA01B: find frequent words
func FrequentWords(pattern string, k int) []string {
	var most_freq_words []string
	dict := make(map[string]int)
	var max_val int
	for i := 0; i < len(pattern)-k+1; i++ {
		frag := pattern[i : i+k]
		if val, ok := dict[frag]; ok {
			new_val := val + 1
			dict[frag] = new_val
			if new_val > max_val {
				max_val = new_val
			}
		} else {
			dict[frag] = 1
		}
	}
	for k, v := range dict {
		if v == max_val {
			most_freq_words = append(most_freq_words, k)
		}
	}
	return most_freq_words
}

// BA01B: find frequent words, fater version
func FreqentWordFaster(pattern string, k int) []string {
	var freq_patterns []string             // prepare return value
	freq_arr := computing_freq(pattern, k) // get frequency array
	// get maximum frequency in frequency array
	max_freq := 0
	for _, v := range freq_arr {
		if v > max_freq {
			max_freq = v
		}
	}
	// extract indices of max frequency, then number -> string
	for i, v := range freq_arr {
		if v == max_freq {
			kmer := num_to_pattern(i, k)
			freq_patterns = append(freq_patterns, kmer)
		}
	}
	return freq_patterns
}

// main
func main() {
	filename := "/home/wsl/rosalind/data/ba01b.txt"
	lines, err := ReadLines(filename)
	if err != nil {
		log.Fatalf("ReadLines: %s", err)
	}

	pattern := lines[0]
	num_k, err := strconv.Atoi(lines[1])
	if err != nil {
		fmt.Println(err)
		os.Exit(2) // 0 means success
	}

	// FrequentWords
	start := time.Now()
	words := FrequentWords(pattern, num_k)
	for _, word := range words {
		fmt.Println(word)
	}
	elapsed := time.Since(start)
	fmt.Printf("page took %s\n", elapsed) // page took 173.9µs

	// FasterFrequentWords
	start = time.Now()
	freq_words := FreqentWordFaster(pattern, num_k)
	for _, v := range freq_words {
		fmt.Println(v)
	}
	elapsed = time.Since(start)
	fmt.Printf("page took %s\n", elapsed) // page took 402.037ms, TODO: Why?
}

// helper fx: readlines
func ReadLines(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	return lines, scanner.Err()
}

// helper fx: print Map
func PrintMap(dict map[string]int) {
	for k, v := range dict {
		fmt.Println(k, ": ", v)
	}
}

// helper fx: integer Pow fx
func int_pow(base int, pw int) int {
	val := math.Pow(float64(base), float64(pw))
	return int(val)
}

// BA01L: pattern to number
func pattern_to_num(pattern string) int {
	var number int
	len := len(pattern)
	for i := 0; i < len; i++ {
		pw := len - i - 1
		var nuc int
		switch pattern[i] {
		case 'A':
			nuc = 0
		case 'C':
			nuc = 1
		case 'G':
			nuc = 2
		case 'T':
			nuc = 3
		}
		pw_val := int_pow(4, pw)
		number = number + (pw_val * nuc)
	}
	return number
}

// BA01M: number to pattern
func num_to_pattern(number int, k int) string {
	var pattern string
	for i := 0; i < k; i++ {
		temp := number / int_pow(4, k-i-1)
		switch temp {
		case 0:
			pattern += "A"
		case 1:
			pattern += "C"
		case 2:
			pattern += "G"
		case 3:
			pattern += "T"
		}
		number -= temp * int_pow(4, k-i-1)
	}
	return pattern
}

// BA01K: compute frequency
func computing_freq(pattern string, k int) []int {
	// create an empty slice of size 4^k
	freq_arr := make([]int, int_pow(4, k))
	for i := range freq_arr {
		freq_arr[i] = 0
	}

	// kmer -> number, update freq_arr
	var kmer string
	var number int
	for i := 0; i < len(pattern)-k+1; i++ {
		kmer = pattern[i : i+k]
		number = pattern_to_num(kmer)
		freq_arr[number] += 1
	}
	return freq_arr
}
