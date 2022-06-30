/*
Rosalind: BA3E (difficulty: 1/5)
Construct the "De Bruijn Graph" of a Collection of k-mers

Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times),
we define "CompositionGraph(Patterns)" as a graph with |Patterns| isolated edges.
Every edge is labeled by a k-mer from Patterns,
and the starting and ending nodes of an edge are labeled
by the prefix and suffix of the k-mer labeling that edge.
We then define the "de Bruijn graph" of Patterns, denoted "DeBruijn(Patterns)",
by gluing identically labeled nodes in CompositionGraph(Patterns),
which yields the following algorithm.

    DEBRUIJN(Patterns)
        represent every k-mer in Patterns as an isolated edge between its prefix and suffix
        glue all nodes with identical labels, yielding the graph DeBruijn(Patterns)
        return DeBruijn(Patterns)

"De Bruijn Graph" from k-mers Problem
Construct the "de Bruijn graph" from a collection of k-mers.

Given: A collection of k-mers Patterns.

Return: The "de Bruijn graph" "DeBruijn(Patterns)", in the form of an adjacency list.

Sample Dataset
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG

Sample Output
AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> PREV: Construct De Bruijn Graph (BA3D)
      -> HERE: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> NEXT: Eulerian Cycle (BA3F)

═════════════════════════════════════════════════

References:
-
*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"time"
)

// BA03E: De Bruijn graph from k-mer patterns
func DeBruijnFromKmers(patterns []string) map[string][]string {
	graph := make(map[string][]string)
	for _, pattern := range patterns {
		prefix := pattern[:len(pattern)-1]
		suffix := pattern[1:]
		if values, ok := graph[prefix]; ok {
			values = append(values, suffix)
			graph[prefix] = values
		} else {
			graph[prefix] = []string{suffix}
		}
	}
	return graph
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03e.txt")

	// print output to screen
	start := time.Now()
	for key, values := range DeBruijnFromKmers(lines) {
		fmt.Print(key, " -> ")
		fmt.Println(strings.Join(values, ","))
	}
	fmt.Println("Execution time: ", time.Since(start))

	/*
		// print output to text file
		old := os.Stdout // save os.Stdout
		f, err := os.Create("/home/wsl/rosalind/data/ba03e_output.txt")
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		os.Stdout = f // redirect os.Stdout to file
		for key, values := range DeBruijnFromKmers(lines) {
			fmt.Print(key, " -> ")
			fmt.Println(strings.Join(values, ","))
		}
		os.Stdout = old // recover os.Stdout
	*/
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
