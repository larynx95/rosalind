/************************************************
Rosalind: BA1L
Implement PatternToNumber

Implement PatternToNumber
Convert a DNA string to a number.

Given: A DNA string Pattern.
Return: PatternToNumber(Pattern).

Sample Dataset
AGT

Sample Output
11

═════════════════════════════════════════════════
Pseudocode:

  ╔════════════════════════════════════════════════════════════════════╗
  ║  PATTERNTONUMBER(Pattern)                                          ║
  ║      if Pattern contains no symbols                                ║
  ║          return 0                                                  ║
  ║      symbol LASTSYMBOL(Pattern)                                    ║
  ║      Prefix PREFIX(Pattern)                                        ║
  ║      return 4 · PATTERNTONUMBER(Prefix) + SYMBOLTONUMBER(symbol)   ║
  ║                                                                    ║
  ╚════════════════════════════════════════════════════════════════════╝

-------------------------------------------------
Plan 1.
- find a rule between given string and index
kmer    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
index    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
x % 4    0  1  2  3  0  1  2  3  0  1  2  3  0  1  2  3
x / 4    0  0  0  0  1  1  1  1  2  2  2  2  3  3  3  3

              A:0, C:1, G:2, T:3

  exp   1 0              2 1 0
  Nuc   T T : 3 3        G A T : 2 0 3
        │ └── 4^0*3      │ │ └── 4^0*2
        └──── 4^1*3      │ └──── 4^1*0
               = 15      └────── 4^2*3
                                  = 50

═════════════════════════════════════════════════

Topic:
- modulo operator
- number types: long long int, long double
- std::pow(), return type (Be careful)
************************************************/

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"time"
)

func PatternToNumberBinary(s string) int64 {
	var o int64
	for _, b := range s {
		o = o<<2 | int64(b&6^b>>1&2)
	}
	return o >> 1
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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01l.txt")
	pattern := lines[0]

	start := time.Now()
	fmt.Println(PatternToNumberBinary(pattern))
	elapsed := time.Since(start)
	fmt.Printf("page took %s\n", elapsed) // 20.5 us

	start = time.Now()
	fmt.Println(PatternToNumber(pattern))
	elapsed = time.Since(start)
	fmt.Printf("page took %s\n", elapsed) // 6.8 us  three times faster

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
