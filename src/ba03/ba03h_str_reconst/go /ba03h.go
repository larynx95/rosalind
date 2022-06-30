/*
Rosalind: BA3H
Reconstruct a String from its k-mer Composition

String Reconstruction Problem
Reconstruct a string from its k-mer composition.

Given: An integer k followed by a list of k-mers Patterns.

Return: A string Text with k-mer composition equal to Patterns.
(If multiple answers exist, you may return any one.)

Sample Dataset
4
CTTA
ACCA
TACC
GGCT
GCTT
TTAC

Sample Output
GGCTTACCA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * PREV: Construct OverapGraph (BA3C)
      -> HERE: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> NEXT: Construct De Bruijn Graph (BA3D)

Plan 1.
  * multiple strings -> overlap graph -> Hamiltonian path
  * split k-mer to two (k-1)-mers (prefix and suffix)
    CTTA -> CTT,TTA
    ACCA -> ACC,CCA
    TACC -> TAC,ACC
    GGCT -> GGC,GCT
    GCTT -> GCT,CTT
    TTAC -> TTA,TAC
    prefixes: [CCT,ACC,TAC,GGC,GCT,TTA]
    suffixes: [TTA,CCA,ACC,GCT,CTT,TAC]

Plan 2.
  * using Eulerian path (BA3G) algorithm
  * multiple strings -> De Bruijn graph -> Eulerian path

═════════════════════════════════════════════════

References:
-
*/

package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"strings"
	"time"
)

func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03h.txt")
	// k, _ := strconv.Atoi(lines[0])
	patterns := lines[1:]
	graph := DeBruijnFromKmers(patterns)
	path := EulerianPath(graph)
	answer := path[0]
	for i := 1; i < len(path); i++ {
		answer += string(path[i][len(path[i])-1])
	}
	println(answer)
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

// helper fx: read graph
func ReadGraph(lines []string) map[string][]string {
	graph := make(map[string][]string)
	for _, line := range lines {
		splitted := strings.Split(line, " -> ")
		key := splitted[0]
		values := strings.Split(splitted[1], ",")
		graph[key] = values
	}
	return graph
}

// helper fx: print graph
func PrintGraph(graph map[string][]string) {
	for key, values := range graph {
		fmt.Print(key, " -> ")
		fmt.Println(strings.Join(values, ","))
	}
}

// print slice
func PrintVec(vec []string) {
	for i, val := range vec {
		if i == len(vec)-1 {
			fmt.Println(val)
		} else {
			fmt.Print(val, "->")
		}
	}
}

// print degrees
func PrintDegrees(degrees map[string][]int) {
	for key, values := range degrees {
		fmt.Print(key, " -> ")
		for i := 0; i < len(values); i++ {
			if i == 1 {
				fmt.Println(values[i])
			} else {
				fmt.Print(values[i], ",")
			}
		}
	}
}

// walk
func Walk(graph map[string][]string, source string) []string {
	rand.Seed(time.Now().UnixNano())
	cycle := []string{source}
	for len(graph[source]) != 0 {
		values := graph[source]
		target := values[rand.Intn(len(values))]
		cycle = append(cycle, target)
		new_values := []string{}
		for _, val := range values {
			if val != target {
				new_values = append(new_values, val)
			}
		}
		graph[source] = new_values
		source = target
	}
	return cycle
}

// check if all values are empty slices
func is_depleted(graph map[string][]string) bool {
	for _, value := range graph {
		if len(value) != 0 {
			return false
		}
	}
	return true
}

// BA03F: Eulerian cycle
func EulerianCycle(graph map[string][]string) []string {
	// select a start node randomly
	rand.Seed(time.Now().UnixNano())
	var nodes []string
	for key := range graph {
		nodes = append(nodes, key)
	}
	start := nodes[rand.Intn(len(nodes))]
	// get a cycle
	cycle := Walk(graph, start)
	// iteration
	for !is_depleted(graph) {
		// select a new start node for the new cycle
		var overlap_idx int
		for idx, node := range cycle {
			if len(graph[node]) != 0 {
				start = node
				overlap_idx = idx
				break
			}
		}
		// get the new cycle
		new_cycle := Walk(graph, start)
		// merge two cycles together (Be careful! Slice operation can return unexpected result.)
		var rest []string
		rest = append(rest, cycle[overlap_idx+1:]...)
		cycle = append(cycle[:overlap_idx], new_cycle...)
		cycle = append(cycle, rest...)
	}
	return cycle
}

// get degrees
func GetDegrees(graph map[string][]string) map[string][]int {
	degrees := make(map[string][]int)
	for key, values := range graph {
		// out-degrees
		if _, ok := degrees[key]; ok { // if key exists in map
			degrees[key][1] += len(values)
		} else { // if key doesn't exit in map
			degrees[key] = []int{0, len(values)}
		}
		// in-degrees
		for _, value := range values {
			if _, ok := degrees[value]; ok { // if key exists in map
				degrees[value][0] += 1
			} else { // if key doesn't exist in map
				degrees[value] = []int{1, 0}
			}
		}
	}
	return degrees
}

// find start, sink node
func FindStartSink(graph map[string][]string) (string, string) {
	var start, sink string
	for key, values := range GetDegrees(graph) {
		if values[0] < values[1] { // if out-degrees > in-degrees
			start = key
		} else if values[0] > values[1] {
			sink = key
		}
	}
	return start, sink
}

// BA03G: Eulerian path
func EulerianPath(graph map[string][]string) []string {
	// connect sink to start node
	start, sink := FindStartSink(graph)
	graph[sink] = append(graph[sink], start)
	// get Eulerian cycle
	cycle := EulerianCycle(graph)
	var idx int
	for i, node := range cycle {
		if node == sink {
			idx = i
		}
	}
	// rearrange cycle: start -> ... -> sink
	if cycle[0] == sink {
		// error case: if sink node is in front of cycle ...
		return cycle[1:]
	}
	path := cycle[idx+1 : len(cycle)-1]
	path = append(path, cycle[:idx+1]...)
	return path
}
