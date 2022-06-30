/*
Rosalind: BA3F (difficulty: 4/5)
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

╔══════════════════════════════════════════════════════════════════════╗
║ EULERIANCYCLE(Graph)                                                 ║
║     form a cycle Cycle by randomly walking in Graph                  ║
║         (don't visit the same edge twice!)                           ║
║     while there are unexplored edges in Graph                        ║
║         select a node newStart in Cycle with still unexplored edges  ║
║         form Cycle' by traversing Cycle (starting at newStart)       ║
║           and then randomly walking                                  ║
║         Cycle <- Cycle'                                              ║
║     return Cycle                                                     ║
╚══════════════════════════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║
║         v <- a node with indegree larger than 1 in G        ║
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║
║                 if NewGraph is connected                    ║
║                     add NewGraph to ALlGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝

Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An "Eulerian directed graph", in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3
1 -> 0
2 -> 1,6
3 -> 2
4 -> 2
5 -> 4
6 -> 5,8
7 -> 9
8 -> 7
9 -> 6

Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

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
      -> PREV: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> HERE: Eulerian Cycle (BA3F)
      -> NEXT: Eulerian Path (BA3G)

Info.
1. Hierholzer's algorithm
  * steps:
    ; just walk through the graph, and get a circle1
    ; if the circle1 is complete whole large cycle without remaining nodes in graph, that's Eulerian cycle
    ; if there're remaining unvisited node in graph yet, just do another walking, get another cycle2
      but in this time, start node should be a node in previous cycle, and it has one or more out-degrees
    ; patch cycle2 to cycle1, repeat previous steps
              4 < 5             4 < 5             4 * 5             4 * 5
              v   ^             v   ^             *   *             *   *
          1 < 2 > 6 > 8     1 *[2]> 6 > 8     1 * 2 * 6 > 8     1 * 2 *[6]* 8
          v   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
          0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
                            10321             10321             103265421
                                                 26542              68796
                                              103265421         1032687965421

═════════════════════════════════════════════════

References:
- Slice length and capacity
  https://go.dev/tour/moretypes/11
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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03f.txt")

	graph := ReadGraph(lines)
	cycle := EulerianCycle(graph)
	fmt.Println(strings.Join(cycle, "->"))
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

// code by others, by sonia
// Graph must be Eulerian, m must be number of edges in Graph.
func EulerianCycle_sonia(Graph [][]int, m int) []int {
	// low end of p is stack of unfinished nodes
	// high end is finished path
	p := make([]int, m+1) // stack + path
	for s := 0; s >= 0; {
		n := p[s]
		if arcs := Graph[n]; len(arcs) > 0 {
			w := arcs[0] // follow first arc
			s++          // push followed node on stack
			p[s] = w
			Graph[n] = arcs[1:] // consume arc
			n = w
		} else { // stuck, move node
			p[m] = n
			s--
			m--
		}
	}
	return p
}
