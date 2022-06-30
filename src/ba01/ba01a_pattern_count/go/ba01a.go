/*
Rosalind: BA1A
Compute the Number of Times a Pattern Appears in a Text

This is the first problem in a collection of "code challenges"
to accompany Bioinformatics Algorithms:
An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

A k-mer is a string of length k.
We define Count(Text, Pattern) as the number of times
that a k-mer Pattern appears as a substring of Text.
For example,

Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.

We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2)
since we should account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern),
our plan is to "slide a window" down Text,
checking whether each k-mer substring of Text matches Pattern.
We will therefore refer to the k-mer starting at position i of Text as Text(i,
k). Throughout this book, we will often use 0-based indexing, meaning that we
count starting at 0 instead of 1. In this case, Text begins at position 0 and
ends at position |Text| − 1 (|Text| denotes the number of symbols in Text). For
example, if Text = GACCATACTG, then Text(4, 3) = ATA. Note that the last k-mer
of Text begins at position |Text| − k, e.g., the last 3-mer of GACCATACTG starts
at position 10 −3 = 7. This discussion results in the following pseudocode for
computing Count(Text, Pattern).

╔═════════════════════════════════════════╗
║ PatternCount(Text, Pattern)             ║
║     count <- 0                          ║
║     for i <- 0 to |Text| - |Pattern|    ║
║         if Text(i, |Pattern|) = Pattern ║
║             count <- count + 1          ║
║     return count                        ║
╚═════════════════════════════════════════╝

Implement PatternCount
Given: {DNA strings}} Text and Pattern.
Return: Count(Text, Pattern).

Sample Dataset
GCGCG
GCG

Sample Output
2

*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"time"
)

// BA01A: PatternCount
func PatternCount(text, pattern string) int {
	count := 0
	for i := 0; i < len(text)-len(pattern)+1; i++ {
		sub := text[i : i+len(pattern)]
		if sub == pattern {
			count++
		}
	}
	return count
}

// main
func main() {
	var filename string = "/home/wsl/rosalind/data/ba01a.txt"
	lines, err := ReadLines(filename)
	if err != nil {
		log.Fatalf("lines: %s", err)
	}
	text := lines[0]
	pattern := lines[1]

	start := time.Now()
	println(PatternCount(text, pattern))
	fmt.Printf("page took %s\n", time.Since(start))
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
