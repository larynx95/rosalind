/*
Rosalind: BA1D
Find All Occurrences of a Pattern in a String

In this problem, we ask a simple question: how many times can one string occur
as a substring of another? Recall from “Find the Most Frequent Words in a
String” that different occurrences of a substring can overlap with each other.
For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem
Find all occurrences of a pattern in a string.

Given: Strings Pattern and Genome.

Return: All starting positions in Genome where Pattern appears as a substring.
Use 0-based indexing.

Sample Dataset
ATAT
GATATATGCATATACTT

Sample Output
1 3 9

═════════════════════════════════════════════════

Topics:
- string comparison
- substring
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
)

// BA01D: find all occurrence
func all_occur_indices(pattern string, genome string) []int {
	var indices []int
	for i := 0; i < len(genome)-len(pattern)+1; i++ {
		if pattern == genome[i:i+len(pattern)] {
			indices = append(indices, i)
		}
	}
	return indices
}

// main
func main() {
	lines, err := read_lines("/home/wsl/rosalind/data/ba01d.txt")
	if err != nil {
		log.Fatalf("read_lines: %s", err)
	}
	pattern := lines[0]
	genome := lines[1]

	indices := all_occur_indices(pattern, genome)
	for _, v := range indices {
		fmt.Print(v, " ")
	}
	fmt.Println()
}

// helper fx: readlines
func read_lines(path string) ([]string, error) {
	// open file, close file
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
