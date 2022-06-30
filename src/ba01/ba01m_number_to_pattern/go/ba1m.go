/************************************************
Rosalind: BA1M
Implement NumberToPattern

Convert an integer to its corresponding DNA string.

Given: Integers index and k.
Return: NumberToPattern(index, k).

Sample Dataset
45
4

Sample Output
AGTC

═════════════════════════════════════════════════
Pseudocode:

  ╔══════════════════════════════════════════════════════════╗
  ║  NUMBERTOPATTERN(index , k)                              ║
  ║      if k = 1                                            ║
  ║          return NUMBERTOSYMBOL(index)                    ║
  ║      prefixIndex QUOTIENT(index, 4)                      ║
  ║      r REMAINDER(index, 4)                               ║
  ║      symbol NUMBERTOSYMBOL(r)                            ║
  ║      PrefixPattern NUMBERTOPATTERN(prefixIndex, k  1)    ║
  ║      return concatenation of PrefixPattern with symbol   ║
  ║                                                          ║
  ╚══════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

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
- '/', '%' operators
- TODO: Implement recursive number_to_pattern fx.
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

// code by sonia
// too difficult for me, and... strangely... slower than my code (Why?)
// TODO: Analyze this code
func NumberToPatternBinary(n int64, k int) string {
	a := make([]byte, k)
	for k > 0 {
		k--
		a[k] = "ACGT"[n&3]
		n >>= 2
	}
	return string(a)
}

// BA01M: NumberToPattern
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

// main
func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba01m.txt")
	number, _ := strconv.Atoi(lines[0])
	k, _ := strconv.Atoi(lines[1])

	// sonia's code
	// Bitwise operation is known to be very fast cuz it is machine-friendly.
	// But this took seven times slower than my elementary (newbie's) code. Why?
	start := time.Now()
	fmt.Println(NumberToPatternBinary(int64(number), k))
	elapsed := time.Since(start)
	fmt.Printf("Execution time (bitwise, machine friendly operation): %s\n", elapsed) // 29.1 us

	// my code
	start = time.Now()
	fmt.Println(NumberToPattern(number, k))
	elapsed = time.Since(start)
	fmt.Printf("Execution time (Pythonic code): %s\n", elapsed) // 3.9 us
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
