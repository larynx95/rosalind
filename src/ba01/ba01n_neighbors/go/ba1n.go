/************************************************
Rosalind: BA1N
Generate the d-Neighborhood of a String

The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming
distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String
Find all the neighbors of a pattern.

Given: A DNA string Pattern and an integer d.
Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
ACG
1

Sample Output
CCG
TCG
GCG
AAG
ATG
AGG
ACA
ACC
ACT
ACG

═════════════════════════════════════════════════

Pseudocode:
  * immediate neighbors
    ACG, d=1
    A --> CCG, GCG, TCG
    C --> AAG, AGG, ATG
    G --> ACA, ACC, ACT
    ╔═══════════════════════════════════════════════════════════════════════╗
    ║  IMMEDIATENEIGHBORS(Pattern)                                          ║
    ║    Neighborhood <- the set consisting of the single string Pattern    ║
    ║    for i = 1 to |Pattern|                                             ║
    ║      symbol <- i-th nucleotide of Pattern                             ║
    ║      for each nucleotide x different from symbol                      ║
    ║        Neighbor <- Pattern with the i-th nucleotide substituted by x  ║
    ║        add Neighbor to Neighborhood                                   ║
    ║    return Neighborhood                                                ║
    ╚═══════════════════════════════════════════════════════════════════════╝

  * iterative neighbors
    ╔══════════════════════════════════════════════════════════════╗
    ║  ITERATIVENEIGHBORS(Pattern, d)                              ║
    ║    Neighborhood <- set consisting of single string Pattern   ║
    ║    for j = 1 to d                                            ║
    ║      for each string Pattern' in Neighborhood                ║
    ║        add IMMEDIATENEIGHBORS(Pattern') to Neighborhood      ║
    ║        remove duplicates from Neighborhood                   ║
    ║    return Neighborhood                                       ║
    ╚══════════════════════════════════════════════════════════════╝

  * neighbors - recursive function
    (1) base case
      - if d == 0          => {pattern}
      - if |pattern| == 0  => {}
      - if |pattern| == 1  => {'A', 'C', 'G', 'T'}
    (2) recursive case
      - if neighbors(pattern[1:]) < d   => can prepend any nucleotide
      - if neighbors(pattern[1:]) >= d  => can change the first nucleotide
    ╔══════════════════════════════════════════════════════════╗
    ║  NEIGHBORS(Pattern, d)                                   ║
    ║    if d = 0                                              ║
    ║      return {Pattern}                                    ║
    ║    if |Pattern| = 1                                      ║
    ║      return {A, C, G, T}                                 ║
    ║    Neighborhood <- an empty set                          ║
    ║    SuffixNeighbors <- NEIGHBORS(SUFFIX(Pattern), d)      ║
    ║    for each string Text from SuffixNeighbors             ║
    ║      if HAMMINGDISTANCE(SUFFIX(Pattern), Text) < d       ║
    ║        for each nucleotide x                             ║
    ║          add x + Text to Neighborhood                    ║
    ║      else                                                ║
    ║        add FIRSTSYMBOL(Pattern) + Text to Neighborhood   ║
    ║    return Neighborhood                                   ║
    ╚══════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- How to replace nth char from a string in Go
  https://stackoverflow.com/questions/37688457/how-to-replace-nth-char-from-a-string-in-go
  ; Strings are immutable.
- How to assign string to bytes array
  https://stackoverflow.com/questions/8032170/how-to-assign-string-to-bytes-array
- How can I change the for loop iterator type?
  https://stackoverflow.com/questions/15030327/how-can-i-change-the-for-loop-iterator-type
- How to remove duplicates strings or int from Slice in Go
  https://stackoverflow.com/questions/66643946/how-to-remove-duplicates-strings-or-int-from-slice-in-go
- Slice of Slices in Golang
  https://www.geeksforgeeks.org/slice-of-slices-in-golang/
- string vs. []string as function argument

************************************************/

package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
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

func ImmediateNeighbors_error(pattern string) []string {
	// TODO: fix this
	var neighborhood []string
	neighborhood = append(neighborhood, pattern)
	for i := 0; i < len(pattern); i++ {
		symbol := pattern[i]
		for j := 0; j < len("ACGT"); j++ {
			if pattern[j] != symbol { // <-- problems in (a) index, (b) comparison, Don't do this kind of mistake.
				neighbor := pattern[:i] + string("ACGT"[j]) + pattern[i+1:]
				neighborhood = append(neighborhood, neighbor)
			}
		}
	}
	return neighborhood
}

func ImmediateNeighbors_fixed(pattern string) []string {
	var neighborhood []string
	neighborhood = append(neighborhood, pattern)
	for i := 0; i < len(pattern); i++ {
		symbol := pattern[i]
		nucleotides := []byte("ACGT")
		for j := 0; j < len("ACGT"); j++ {
			nuc := nucleotides[j]
			if symbol != nuc {
				neighbor := pattern[:i] + string(nuc) + pattern[i+1:]
				neighborhood = append(neighborhood, neighbor)
			}
		}
	}
	return neighborhood
}

func ImmediateNeighbors(pattern string) []string {
	var neighborhood []string // neighborhood = []string{pattern}
	neighborhood = append(neighborhood, pattern)
	for i, symbol := range pattern {
		for _, nuc := range "ACGT" {
			if symbol != nuc {
				neighbor := pattern[:i] + string(nuc) + pattern[i+1:]
				neighborhood = append(neighborhood, neighbor)
			}
		}
	}
	return neighborhood
}

// BA01N: neighbors
func Neighbors_iter(pattern string, distance int) []string {
	var neighborhood []string
	neighborhood = append(neighborhood, pattern)
	for i := 0; i < distance; i++ {
		for _, pat := range neighborhood { // 'neighborhood' changes every loop.
			neighborhood = append(neighborhood, ImmediateNeighbors(pat)...)
			neighborhood = RemoveDuplicateStr(neighborhood) // no 'Set' in Golnag
		}
	}
	sort.Strings(neighborhood) // package 'sort.Strings([]string)'
	return neighborhood
}

// BA01N: neighbors
func Neighbors_rec(pattern string, distance int) []string {
	if distance == 0 {
		return []string{pattern}
	}
	if len(pattern) == 1 {
		return []string{"A", "C", "G", "T"}
	}
	neighborhood := []string{}
	suffixPattern := pattern[1:]
	suffixNeighbors := Neighbors_rec(suffixPattern, distance)
	for _, neighbor := range suffixNeighbors {
		if HammingDistance(suffixPattern, neighbor) < distance {
			for _, nuc := range []string{"A", "C", "G", "T"} {
				neighborhood = append(neighborhood, string(nuc)+neighbor)
			}
		} else {
			neighborhood = append(neighborhood, string(pattern[0])+neighbor)
		}
	}
	sort.Strings(neighborhood)
	return neighborhood
}

// code by sonia
// too difficult for me
func Neighbors1(s string, h int) []string {
	v := []string{s}        // result list
	const sym = "A C T G"   // used for encoding symbols from indexes
	buf := []byte(s)        // buffer for generating neighbors
	var f func([]byte, int) // recursive function
	// closure
	f = func(t []byte, h int) {
		for i := 0; i < len(t); i++ {
			sub := t[i:]
			b := sub[0]
			for vb := byte(0); ; vb += 2 {
				if vb == b&6 {
					vb += 2
				}
				if vb == 8 {
					break
				}
				sub[0] = sym[vb]
				v = append(v, string(buf))
				if h > 1 && len(sub) > 1 {
					f(sub[1:], h-1)
				}
			}
			sub[0] = b
		}
	}
	f(buf, h)
	return v
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01n.txt")
	pattern := lines[0]
	k, _ := strconv.Atoi(lines[1])

	// check execution time, ImmediateNeighbors fx
	start := time.Now()
	nbs_imm := ImmediateNeighbors_fixed(pattern)
	for _, kmer := range nbs_imm {
		fmt.Println(kmer)
	}
	fmt.Printf("count: %d\n", len(nbs_imm))
	elapsed := time.Since(start)
	fmt.Printf("Execution time (Immediate Neighbors): %s\n\n", elapsed)

	// Negibors1 fx
	start = time.Now()
	nbs_1 := Neighbors1(pattern, k)
	for _, kmer := range nbs_1 {
		fmt.Println(kmer)
	}
	fmt.Printf("count: %d\n", len(nbs_1))
	elapsed = time.Since(start)
	fmt.Printf("Execution time (Neighbors, closure): %s\n\n", elapsed)

	// Neighbors_iter fx.
	start = time.Now()
	nbs_iter := Neighbors_iter(pattern, 2)
	for _, kmer := range nbs_iter {
		fmt.Println(kmer)
	}
	fmt.Printf("count: %d\n", len(nbs_iter))
	elapsed = time.Since(start)
	fmt.Printf("Execution time (Neighbors, iterative): %s\n\n", elapsed)

	// Neighbors_rec fx.
	start = time.Now()
	nbs_rec := Neighbors_rec(pattern, 2)
	for _, kmer := range nbs_rec {
		fmt.Println(kmer)
	}
	fmt.Printf("count: %d\n", len(nbs_rec))
	elapsed = time.Since(start)
	fmt.Printf("Execution time (Neighbors, recursive): %s\n\n", elapsed)

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

// helper fx: remove duplicates (string type)
// https://stackoverflow.com/questions/66643946/how-to-remove-duplicates-strings-or-int-from-slice-in-go
// get unique keys from map
func RemoveDuplicateStr(strSlice []string) []string {
	allKeys := make(map[string]bool)
	list := []string{}
	for _, item := range strSlice {
		if _, value := allKeys[item]; !value {
			allKeys[item] = true
			list = append(list, item)
		}
	}
	return list
}

// helper fx: remove duplicates (integer type)
// https://stackoverflow.com/questions/66643946/how-to-remove-duplicates-strings-or-int-from-slice-in-go
func removeDuplicateInt(intSlice []int) []int {
	allKeys := make(map[int]bool)
	list := []int{}
	for _, item := range intSlice {
		if _, value := allKeys[item]; !value {
			allKeys[item] = true
			list = append(list, item)
		}
	}
	return list
}

// helper fx: another remove duplicates fx
// sort string clice, then remove repeating elements
func RemoveDuplicateStrSorting(s []string) []string {
	if len(s) < 1 {
		return s
	}
	sort.Strings(s)
	prev := 1
	for curr := 1; curr < len(s); curr++ {
		if s[curr-1] != s[curr] {
			s[prev] = s[curr]
			prev++
		}
	}
	return s[:prev]
}

// generic
type SliceType interface {
	~string | ~int | ~float64 // add more *comparable* types as needed
}

func removeDuplicates[T SliceType](s []T) []T {
	if len(s) < 1 {
		return s
	}
	// sort
	sort.SliceStable(s, func(i, j int) bool {
		return s[i] < s[j]
	})
	prev := 1
	for curr := 1; curr < len(s); curr++ {
		if s[curr-1] != s[curr] {
			s[prev] = s[curr]
			prev++
		}
	}
	return s[:prev]
}
