/************************************************
Rosalind: BA1K
Generate the Frequency Array of a String

Given an integer k, we define the frequency array of a string Text as an array
of length 4k, where the i-th element of the array holds the number of times that
the i-th k-mer (in the lexicographic order) appears in Text (see Figure 1.

kmer      AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
index      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
frequency  3   0   2   0   1   0   0   0   0   1   3   1   0   0   1   0

Computing a Frequency Array
Generate the frequency array of a DNA string.

Given: A DNA string Text and an integer k.

Return: The frequency array of k-mers in Text.

Sample Dataset
ACGCGGCTCTGAAA
2

Sample Output
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

═════════════════════════════════════════════════
Pseudocode:

  ╔═══════════════════════════════════════════════════╗
  ║  COMPUTINGFREQUENCIES(Text, k)                    ║
  ║    for i <- 0 to 4^k- 1                           ║
  ║      FREQUENCYARRAY(i) <- 0                       ║
  ║    for i <- 0 to |Text| - k                       ║
  ║      Pattern <- Text(i, k)                        ║
  ║      j <- PATTERNTONUMBER(Pattern)                ║
  ║      FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1   ║
  ║    return FREQUENCYARRAY                          ║
  ╚═══════════════════════════════════════════════════╝

═════════════════════════════════════════════════
References:
- Go << and >> operators
  https://stackoverflow.com/questions/5801008/go-and-operators
  n << x : n * 2^x
  y >> z : y / 2^z
- How to check if a map contains a key in Go?
  https://stackoverflow.com/questions/2050391/how-to-check-if-a-map-contains-a-key-in-go
  if val, ok := dict["foo"]; ok {
    //do something here
  }
- What is the difference between "1 << x" and "pow(2, x)"?
  https://stackoverflow.com/questions/47416967/what-is-the-difference-between-1-x-and-pow2-x
- Iterating through a golang map
  https://stackoverflow.com/questions/8018719/iterating-through-a-golang-map

************************************************/

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"
)

// BA01K: ComputeFrequency, as integer slice
func ComputeFrequency(text string, k int) []int {
	var freqArr = make([]int, int(math.Pow(4, float64(k)))) // empty (all zero) slice
	for i := 0; i < len(text)-k+1; i++ {
		num := PatternToNumber(text[i : i+k])
		freqArr[num] += 1
	}
	return freqArr
}

// BA01K: ComputeFrequency, as Map
func ComputeFrequencyMap(text string, k int) map[int]int {
	var freqMap = make(map[int]int)
	for i := 0; i < len(text); i++ {
		num := PatternToNumber(text[i : i+k])
		if _, ok := freqMap[num]; ok {
			freqMap[num] += 1
		} else {
			freqMap[num] = 1
		}
	}
	return freqMap
}

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01k.txt")
	text := lines[0]
	k, _ := strconv.Atoi(lines[1])

	// ComputeFrequency by butwise operation
	start := time.Now()
	for _, f := range ComputeFrequencyBitwise(text, k) {
		fmt.Print(f, " ")
	}
	fmt.Println()
	elapsed := time.Since(start)
	fmt.Printf("Execution time (bitwise): %s\n\n", elapsed)

	// ComputeFrequency by algorithm in textbook
	start = time.Now()
	for _, f := range ComputeFrequency(text, k) {
		fmt.Print(f, " ")
	}
	fmt.Println()
	elapsed = time.Since(start)
	fmt.Printf("Execution time (textbook algorithm): %s\n\n", elapsed)

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

// BA01L: pattern to number
func PatternToNumber(pattern string) int {
	var num = 0
	var l = len(pattern)
	for i, nuc := range pattern {
		temp := int(math.Pow(4, float64(l-i-1)))
		if nuc == 'A' {
			num += temp * 0
		} else if nuc == 'C' {
			num += temp * 1
		} else if nuc == 'G' {
			num += temp * 2
		} else {
			num += temp * 3
		}
	}
	return num
}

// BA01M: number to pattern
func NumberToPattern(number int, k int) string {
	var pattern []byte
	for i := 0; i < k; i++ {
		temp := number / int(math.Pow(4, float64(k-i-1)))
		if temp == 0 {
			pattern = append(pattern, 'A')
		} else if temp == 1 {
			pattern = append(pattern, 'C')
		} else if temp == 2 {
			pattern = append(pattern, 'G')
		} else if temp == 3 {
			pattern = append(pattern, 'T')
		}
		number -= temp * int(math.Pow(4, float64(k-i-1)))
	}
	return string(pattern)
}

/* code by sonia - analyze this code later */

// TODO: Analyze this code!
func PatternToNumberBitwise(s string) int64 {
	var o int64
	for _, b := range s {
		o = o<<2 | int64(b&6^b>>1&2) // bitwise operation OR, XOR
	}
	return o >> 1
}

func ComputeFrequencyBitwise(Text string, k int) []int {
	a := make([]int, 1<<(2*uint(k))) // bitwise operation 2*2^k == 4^k
	n := PatternToNumberBitwise(Text[:k])
	a[n] = 1
	mask := int64(len(a) - 1)
	for k < len(Text) {
		b := Text[k]
		n = n<<2&mask | int64(b>>1&3^b>>2&1)
		a[n]++
		k++
	}
	return a
}
