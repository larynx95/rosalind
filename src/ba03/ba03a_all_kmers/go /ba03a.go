/*
 */

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"time"
)

// BA03A: get all k-mer patterns
func AllKmers(text string, k int) []string {
	kmers := []string{}
	for i := 0; i < len(text)-k+1; i++ {
		kmers = append(kmers, text[i:i+k])
	}
	return kmers
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03a.txt")
	k, _ := strconv.Atoi(lines[0])
	text := lines[1]

	start := time.Now()
	for _, kmer := range AllKmers(text, k) {
		fmt.Println(kmer)
	}
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
