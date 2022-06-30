/*
Rosalind: BA3D
Construct the "De Bruijn[broin]" Graph of a String

Given a genome Text, PathGraph_k(Text) is the path consisting of |Text| - k + 1 edges,
where the i-th edge of this path is labeled by the i-th k-mer in Text
and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text.
The de Bruijn graph DeBruijn_k(Text) is formed
by gluing identically labeled nodes in PathGraph_k(Text).

De Bruijn Graph from a String Problem
Construct the de Bruijn graph of a string.

Given: An integer k and a string Text.

Return:DeBruijn_k(Text), in the form of an adjacency list.

Sample Dataset
4
AAGATTCTCTAC

Sample Output (1:many)
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC
TTC -> TCT

═════════════════════════════════════════════════

References:
- How to check if a map contains a key in Go?
  https://stackoverflow.com/questions/2050391/how-to-check-if-a-map-contains-a-key-in-go
- golang why don't we have a set datastructure [closed]
  https://stackoverflow.com/questions/34018908/golang-why-dont-we-have-a-set-datastructure
*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"
)

// BA03D: De Bruijn graph from a text
func DeBruijnFromText(text string, k int) map[string][]string {
	graph := make(map[string][]string)
	for i := 0; i < len(text)-k+1; i++ {
		prefix := text[i : i+k-1]
		suffix := text[i+1 : i+k]
		if values, ok := graph[prefix]; ok {
			values = append(values, suffix) // TODO:  can be duplicated
			graph[prefix] = values
		} else {
			graph[prefix] = []string{suffix}
		}
	}
	return graph
}

// BA03D: De Bruijn graph from a text, set version
func DeBruijnFromText_set(text string, k int) map[string](map[string]int) {
	graph := make(map[string]map[string]int)
	for i := 0; i < len(text)-k+1; i++ {
		prefix := text[i : i+k-1]
		suffix := text[i+1 : i+k]
		if _, ok := graph[prefix]; ok {
			graph[prefix][suffix] = 1
		} else {
			graph[prefix] = map[string]int{suffix: 1}
		}
	}
	return graph
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03d.txt")
	k, _ := strconv.Atoi(lines[0])
	//text := lines[1]

	// check time, print output to screen
	start := time.Now()
	for key, values := range DeBruijnFromText("AAGATTCTCTACAAGATTCTCTAC", k) {
		fmt.Print(key, " -> ")
		fmt.Print(strings.Join(values, ","))
		fmt.Println()
	}
	fmt.Println("Execution time: ", time.Since(start))
	println()

	/*
		// print output to text file
		old := os.Stdout // save os.Stdout
		f, err := os.Create("/home/wsl/rosalind/data/ba03d_output.txt")
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		os.Stdout = f // redirect os.Stdout to file
		for source, targets := range DeBruijnFromText_set(text, k) {
			fmt.Print(source, " -> ")
			temp := []string{}
			for target := range targets {
				temp = append(temp, target)
			}
			fmt.Print(strings.Join(temp, ","))
			fmt.Println()
		}
		os.Stdout = old // recover os.Stdout
	*/

	// set version (no duplication)
	for source, targets := range DeBruijnFromText_set("AAGATTCTCTACAAGATTCTCTAC", k) {
		fmt.Print(source, " -> ")
		temp := []string{}
		for target := range targets {
			temp = append(temp, target)
		}
		fmt.Print(strings.Join(temp, ","))
		fmt.Println()
	}
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
