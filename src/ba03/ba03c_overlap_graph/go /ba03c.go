/*
Rosalind: BA3C
Construct the Overlap Graph of a Collection of k-mers

In this chapter, we use the terms prefix and suffix
to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.

Given an arbitrary collection of k-mers Patterns,
we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern'
by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern').
The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

Overlap Graph Problem
Construct the overlap graph of a collection of k-mers.

Given:
A collection Patterns of k-mers.

Return:
The overlap graph Overlap(Patterns), in the form of an adjacency list.

Sample Dataset
ATGCG
GCATG
CATGC
AGGCA
GGCAT

Sample Output (1:1)
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * HERE: Construct OverapGraph (BA3C)
      -> NEXT: Reconstruct a String from its k-mer Composition (BA3H)

Plan 1.
  - deep copy patterns, create two list
  - use nested loop to make adjacent list

═════════════════════════════════════════════════
*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"time"
)

// BA03C: overlap graph
func overlap_graph(patterns []string) map[string][]string {
	dic := make(map[string][]string)
	for _, source := range patterns {
		for _, target := range patterns {
			if source[1:] == target[0:len(target)-1] {
				if val, ok := dic[source]; ok {
					val = append(val, target)
				} else {
					dic[source] = []string{target}
				}
			}
		}
	}
	return dic
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03c.txt")

	start := time.Now()
	for source, targets := range overlap_graph(lines) {
		fmt.Print(source, " -> ")
		for _, target := range targets {
			fmt.Print(target, " ")
		}
		println()
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
