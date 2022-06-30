/*
Rosalind: BA3B
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern[1], ... , Pattern[n]
such that the last k - 1 symbols of Pattern[i] are equal
to the first k - 1 symbols of Pattern[i+1] for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> PREV: Generate the k-mer Composition of a String (BA3A)
      -> HERE: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * NEXT: Construct OverapGraph (BA3C)
*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"time"
)

// BA03B: string spelled by genome path
func str_spelled_by_genome_path(patterns []string) string {
	result := patterns[0]
	for _, s := range patterns[1:] {
		result += string(s[len(s)-1])
	}
	return result
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03b.txt")

	start := time.Now()
	fmt.Println(str_spelled_by_genome_path(lines))
	fmt.Println("Execution time:", time.Since(start))
}

// helper fx: read lines
func ReadLines(path string) []string {
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
