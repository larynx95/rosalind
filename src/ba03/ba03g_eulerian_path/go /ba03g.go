/*
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node) is called an Eulerian path.

Eulerian Path Problem
Find an Eulerian path in a graph.

Given: A directed graph that contains an Eulerian path,
       where the graph is given in the form of an adjacency list.

Return: An Eulerian path in this graph.

Sample Dataset
0 -> 2
1 -> 3
2 -> 1
3 -> 0,4
6 -> 3,7
7 -> 8
8 -> 9
9 -> 6

Sample Output
6->7->8->9->6->3->0->2->1->3->4

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
      -> Construct De Bruijn Graph (BA3D)
      -> Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> PREV: Eulerian Cycle (BA3F)
        -> NEXT: k-universal string problem (BA3I)
      -> HERE: Eulerian Path (BA3G)

Info.
  * writing down graph to figure out what this exercise is:
    {0:[2], 1:[3], 2:[1], 3:[0, 4], 6:[3,7], 7:[8], 8:[9], 9:[6]}
                  4
                  ↑
          9 → 6 → 3 → 0       Not circular!!!
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
  * problems:
    - This graph is not circular graph.
      Can any node be a starting point (source) or an ending point (target)?
    - How can I find the end node?
      a. comparing keys and values in dictionary (De Bruijn graph)
         dict.keys:   0123  6789
         dict.values: 012334 789   --> 4 is the end!
      b. difference bw in-degree and out-degree
         In-degree must be larger than out-degree in the last node.
         0:0, 1:0, 2:0, 3:0, 4:1, 6:-1, 7:0, 8:0, 9:0  --> The node 4 is the last!

Plan 1.
  * re-using "EulerianCycle" function
    - Adding an extra edge from w (end node) to v (start node)
      transforms the Eulerian path into an Eulerian cycle.
    - convert Eulerian path to Eulerian cycle, then re-convert cycle to path

  * out-degree +1, in-degree -1
                  4
                  ↑
          9 → 6 → 3 → 0
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
    unique nodes: 0    1    2    3    4    6    7    8    9
    --------------------------------------------------------
      0:[2]       +1       -1
      1:[3]           +1        -1
      2:[1]           -1   +1
      3:[0,4]     -1            +2   -1
      6:[3,7]                   -1        +2   -1
      7:[8]                                    +1   -1
      8:[9]                                         +1   -1
      9:[6]                               -1             +1
    --------------------------------------------------------
                  0    0    0    0   -1   +1    0    0    0   --> start:6 end (sink):4
  * add edge (end, start) to graph, get Eulerian cycle: 4 -> 6
  * trim the last repetitive elem in cycle, then rotate cycle to get Eulerian path

═════════════════════════════════════════════════

References:
- Concatenate two slices in Go
  https://stackoverflow.com/questions/16248241/concatenate-two-slices-in-go
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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03g.txt")
	graph := ReadGraph(lines)
	path := EulerianPath(graph)
	PrintVec(path)
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
	for key, _ := range graph {
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
