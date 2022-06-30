/************************************************
Rosalind: BA1F
Find a Position in a Genome Minimizing the Skew

Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference
between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi
(Genome) denote the prefix (i.e., initial substring) of Genome of length i. For
example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i
(from 0 to |Genome|).

Sample Dataset
CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

Sample Output
53 97

-------------------------------------------------

Plan 1.
- get GC skew table, find min value
- Find the indices of the values ​​that match the min value

Topics:
- GC skewness
************************************************/

package main

import (
	"bufio"
	"fmt"
	"os"
)

// BA01F: minimum skewness
func MinSkew(genome string) []int {
	var indices []int   // slice for result
	var skew = []int{0} // slice for skewness
	var minskew = 99999 // minimum value of skewness
	for idx, nuc := range genome {
		last := skew[len(skew)-1]
		if nuc == 'G' {
			skew = append(skew, last+1)
		} else if nuc == 'C' {
			skew = append(skew, last-1)
			if last-1 < minskew {
				minskew = last - 1
				indices = []int{idx + 1}
			} else if last-1 == minskew {
				indices = append(indices, idx+1)
			}
		}
	}
	return indices

}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01f.txt")
	genome := lines[0]
	indices := MinSkew(genome)
	for _, val := range indices {
		fmt.Print(val, " ")
	}
	fmt.Println()
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
