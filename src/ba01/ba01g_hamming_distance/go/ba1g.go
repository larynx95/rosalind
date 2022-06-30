/*
Rosalind: BA1G
Compute the Hamming Distance Between Two Strings

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
For example, CGAAT and CGGAC have two mismatches. The number of mismatches
between strings p and q is called the Hamming distance between these strings and
is denoted HammingDistance(p, q).

Hamming Distance Problem
Compute the Hamming distance between two DNA strings.

Given: Two DNA strings.

Return: An integer value representing the Hamming distance.

Sample Dataset
GGGCCGTTGGT
GGACCGTTGAC

Sample Output
3

*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"time"
)

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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01g.txt")
	pat1 := lines[0]
	pat2 := lines[1]

	start := time.Now()
	fmt.Println(HammingDistance(pat1, pat2))
	fmt.Println("Execution time:", time.Since(start))
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
