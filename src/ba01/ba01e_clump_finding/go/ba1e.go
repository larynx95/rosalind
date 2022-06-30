/*
Rosalind: BA1E
Find Patterns Forming Clumps in a String

TODO: Find more efficient algorithms!

Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger)
string Genome if there is an interval of Genome of length L in which Pattern
appears at least t times. For example, TGCA forms a (25,3)-clump in the
following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem
Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset
CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4

Sample Output
CGACA GAAGA AATGT

═════════════════════════════════════════════════
Plan 1.
  - use Map, instead of Array (or Slice) for better performance
  - The algorithm of ClumpFindingBetter fx is very great. Analyze it thouroughly.
  - Golang doesn't have built-in Set datastructure. Use Map instead.

Pseudocode:
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║  CLUMPFINDING(Genome, k, t, L)                                        ║
  ║    FrequentPatterns <- an empty set                                   ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      CLUMP(i) <- 0                                                    ║
  ║    for i <- 0 to |Genome| - L                                         ║
  ║      Text <- the string of length L starting at position i in Genome  ║
  ║      FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)                  ║
  ║      for index <- 0 to 4k - 1                                         ║
  ║        if FREQUENCYARRAY(index) >= t                                  ║
  ║          CLUMP(index) <- 1                                            ║
  ║    for i <- 0 to 4k - 1                                               ║
  ║      if CLUMP(i) = 1                                                  ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                               ║
  ║        add Pattern to the set FrequentPatterns                        ║
  ║    return FrequentPatterns                                            ║
  ╚═══════════════════════════════════════════════════════════════════════╝

  ╔═══════════════════════════════════════════════════════════╗
  ║  CLUMPFINDINGBetter(Genome, k, t, L)                      ║
  ║    FrequentPatterns <- an empty set                       ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      CLUMP(i) <- 0                                        ║
  ║    Text <- Genome(0, L)                                   ║
  ║    FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k)        ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if FREQUENCYARRAY(i) >= t                            ║
  ║        CLUMP(i) <- 1                                      ║
  ║    for i <- 1 to |Genome| - L                             ║
  ║      FirstPattern <- Genome(i - 1, k)                     ║
  ║      index <- PATTERNTONUMBER(FirstPattern)               ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) - 1   ║
  ║      LastPattern <- Genome(i + L - k, k)                  ║
  ║      index <- PATTERNTONUMBER(LastPattern)                ║
  ║      FREQUENCYARRAY(index) <- FREQUENCYARRAY(index) + 1   ║
  ║      if FREQUENCYARRAY(index) >= t                        ║
  ║        CLUMP(index) <- 1                                  ║
  ║    for i <- 0 to 4k - 1                                   ║
  ║      if CLUMP(i) = 1                                      ║
  ║        Pattern <- NUMBERTOPATTERN(i, k)                   ║
  ║        add Pattern to the set FrequentPatterns            ║
  ║    return FrequentPatterns                                ║
  ╚═══════════════════════════════════════════════════════════╝

1. find the generalized rule of frequency array

  "AAGCAAAGGTGGGC"  len=14
   vv-----------                                    first ... last
  "AAGCAAAGGTGGG"   len=13, starting from index 0   AA    ...  -
   "AGCAAAGGTGGGC"  len=13, starting from index 1   -     ...  GC
    -----------^^
    common part!

  (1) compute_freq("AAGCAAAGGTGGG", 2)
                    ^^  <--- minus 1
  (2) compute_freq( "AGCAAAGGTGGGC", 2)
                                ^^  <--- plus 1

      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
      3  0  2  0  1  0  0  0  0  1  3  1  0  0  1  0   --- (1)
     -2  0  2  0  1  0  0  0  0  2+ 3  1  0  0  1  0   --- (2)

2. If we know a frequency array of the first clump,
   we can get the frequency array of the whole genome.

    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT AAGCAAAGGTGGGC   fst  lst
    3  0  2  0  1  0  0  0  0  1  1  1  0  0  1  0  AAGCAAAGGTG
   -2  0  2  0  1  0  0  0  0  1  2+ 1  0  0  1  0   AGCAAAGGTGG     AA   GG
    2  0 -1  0  1  0  0  0  0  1  3+ 1  0  0  1  0    GCAAAGGTGGG    AG   GG
    2  0  1  0  1  0  0  0  0 -1+ 3  1  0  0  1  0     CAAAGGTGGGC   GC   GC
    check the frequency before subtraction!

═════════════════════════════════════════════════
References:
- golang why don't we have a set datastructure [closed]
  https://stackoverflow.com/questions/34018908/golang-why-dont-we-have-a-set-datastructure

*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

// BA01E: ClumpFinding
func ClumpFinding(genome string, k, l, t int) []string {
	var freqPatterns []string
	clump := make(map[int]int)
	for i := 0; i < len(genome)-l+1; i++ {
		freqMap := ComputeFrequencyMap(genome[i:i+l], k)
		for num, freq := range freqMap {
			if freq >= t {
				clump[num] = 1
			}
		}
	}
	for num := range clump {
		freqPatterns = append(freqPatterns, NumberToPattern(num, k))
	}
	return freqPatterns
}

// BA01E: ClumpFindingBetter
func ClumpFindingBetter(genome string, k, l, t int) []string {
	var freqPatterns []string  // result array
	clump := make(map[int]int) // map, instead of array (slice)
	// get the frequency map of the first clump, and create a clump frequency map
	fstText := genome[:l]
	freqMap := ComputeFrequencyMap(fstText, k)
	for i, val := range freqMap {
		if val >= t {
			clump[i] = 1
		}
	}
	// algorithm in textbook
	for i := 1; i < len(genome)-l+1; i++ {
		fstPattern := genome[i-1 : i-1+k]
		idx := PatternToNumber(fstPattern)
		freqMap[idx] = freqMap[idx] - 1
		lstPattern := genome[i+l-k : i+l]
		idx = PatternToNumber(lstPattern)
		freqMap[idx] = freqMap[idx] + 1
		if freqMap[idx] >= t {
			clump[idx] = 1
		}
	}
	for num, val := range clump {
		if val == 1 {
			pattern := NumberToPattern(num, k)
			freqPatterns = append(freqPatterns, pattern)
		}
	}
	return freqPatterns
}

// BA01E: ClumpFinding - code by sonia
func Clumps(Genome string, k, L, t int) map[string]struct{} {
	cs := map[string]struct{}{} // clump set. found clumps.
	pm := map[string][]int{}    // position map. start positions by kmer.
	w := L - k                  // window for start positions
	for i, j := 0, k; j <= len(Genome); i, j = i+1, j+1 {
		kmer := Genome[i:j]
		sp := append(pm[kmer], i)
		pm[kmer] = sp
		if len(sp) >= t && i-sp[len(sp)-t] <= w {
			cs[kmer] = struct{}{}
		}
	}
	return cs
}

// main
func main() {
	// read file, extract data
	lines, err := ReadLines("/home/wsl/rosalind/data/ba01e.txt")
	if err != nil {
		log.Fatalf("read_lines: %s", err)
	}
	genome := lines[0]
	var int_vals = []int{}
	for _, val := range strings.Split(lines[1], " ") {
		j, err := strconv.Atoi(val)
		if err != nil {
			panic(err)
		}
		int_vals = append(int_vals, j)
	}
	k := int_vals[0]
	l := int_vals[1]
	t := int_vals[2]

	// check time: ClumFinding code by other person
	start := time.Now()
	for k := range Clumps(genome, k, l, t) {
		fmt.Print(k, " ")
	}
	fmt.Println()
	fmt.Println("Execution time (code by sonia):", time.Since(start)) // 2.971 ms (fastest)

	// check time: ClumpFinding
	start = time.Now()
	for _, pat := range ClumpFinding(genome, k, l, t) {
		fmt.Print(pat, " ")
	}
	fmt.Println()
	fmt.Println("Execution time (using freqMap):", time.Since(start)) // 2.30725 s

	// check time: ClumpFindingBetter
	start = time.Now()
	for _, pat := range ClumpFindingBetter(genome, k, l, t) {
		fmt.Print(pat, " ")
	}
	fmt.Println()
	fmt.Println("Execution time (better):", time.Since(start)) // 8.3515 ms
}

// helper fx: readlines
func ReadLines(path string) ([]string, error) {
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

// BA01H: HammingDistance
func HammingDistance(pattern1, pattern2 string) int {
	var dist = 0
	for i, _ := range pattern1 {
		if pattern1[i] != pattern2[i] {
			dist += 1
		}
	}
	return dist
}

// BA01L: PatternToNumber
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

// BA01K: ComputeFrequency
func ComputeFrequency(text string, k int) []int {
	var freqArr = make([]int, int(math.Pow(4, float64(k))))
	for i := 0; i < len(text)-k+1; i++ {
		num := PatternToNumber(text[i : i+k])
		freqArr[num] += 1
	}
	return freqArr
}

// BA01K: ComputeFrequency, as Map
func ComputeFrequencyMap(text string, k int) map[int]int {
	var freqMap = make(map[int]int)
	for i := 0; i < len(text)-k+1; i++ {
		num := PatternToNumber(text[i : i+k])
		freqMap[num] += 1 // no need to check a key is in a map
	}
	return freqMap
}
