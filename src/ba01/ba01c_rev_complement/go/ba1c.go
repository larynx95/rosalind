/*
Rosalind: BA1C
Find the Reverse Complement of a String

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C'
and 'G'. Given a nucleotide p, we denote its complementary nucleotide as p. The
reverse complement of a DNA string Pattern = p1…pn is the string rPattern = pn …
p1 formed by taking the complement of each nucleotide in Pattern, then reversing
the resulting string.

For example, the reverse complement of Pattern = "GTCA" is rPattern = "TGAC".

Reverse Complement Problem
Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: rPattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT

═════════════════════════════════════════════════
References:
- How to reverse a string in Go?
  https://stackoverflow.com/questions/1752414/how-to-reverse-a-string-in-go

*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
)

// helper fx: reverse string
func ReverseString(pattern string) string {
	var reversed string
	for i := len(pattern) - 1; i >= 0; i-- {
		reversed = reversed + string(pattern[i])
	}
	return reversed
}

// BA01C: ReverseComplement
func ReverseComplement(pattern string) string {
	var complement_string string
	reversed := ReverseString(pattern)
	for i := 0; i < len(reversed); i++ {
		var nucleotide byte
		switch reversed[i] {
		case 'A':
			nucleotide = 'T'
		case 'C':
			nucleotide = 'G'
		case 'G':
			nucleotide = 'C'
		case 'T':
			nucleotide = 'A'
		}
		complement_string += string(nucleotide)
	}
	return complement_string
}

// main
func main() {
	filename := "/home/wsl/rosalind/data/ba01c.txt"
	lines, err := read_lines(filename)
	if err != nil {
		log.Fatalf("read_lines: %s", err)
	}
	pattern := lines[0]
	fmt.Println(ReverseComplement(pattern))
}

// helper fx: ReadLines
func read_lines(path string) ([]string, error) {
	// open file, close file
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	return lines, scanner.Err()
}

// sample code (by sonia)
// using binary operators --> difficult
func RevC(Pattern string) string {
	rc := make([]byte, len(Pattern))
	rcx := len(rc)
	for _, b := range []byte(Pattern) {
		rcx--
		rc[rcx] = ^b&2>>1*17 | 4 ^ b // <--- binary operations
	}
	return string(rc)
}
