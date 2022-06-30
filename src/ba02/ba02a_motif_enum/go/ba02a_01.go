/************************************************
Rosalind: BA2A
Implement MotifEnumeration

Given a collection of strings Dna and an integer d,
a k-mer is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
The following algorithm finds (k,d)-motifs.

╔════════════════════════════════════════════════════════════════════════════════════╗
║ MOTIFENUMERATION(Dna, k, d)                                                        ║
║     Patterns <- an empty set                                                       ║
║     for each k-mer Pattern in Dna                                                  ║
║         for each k-mer Pattern' differing from Pattern by at most d mismatches     ║
║             if Pattern' appears in each string from Dna with at most d mismatches  ║
║                 add Pattern' to Patterns                                           ║
║     remove duplicates from Patterns                                                ║
║     return Patterns                                                                ║
╚════════════════════════════════════════════════════════════════════════════════════╝

Implanted Motif Problem
Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.

Given: Integers k and d, followed by a collection of strings Dna.

Return: All (k, d)-motifs in Dna.

Sample Dataset
3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output
ATA ATT GTT TTT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Which DNA patterns play the role of molecular clock?
    * find motifs from each of DNA string (e.g. ten 15-mers form ten DNA strings)- different from BA1B, BA1I
      -> "Motif Finding Problem":
         exploring all Motifs n Dna => then deriving the consensus string from Motifs
      ↓
    * Motivation
      -> HERE: finding all (k, d)-motifs appearing in every DNA string (BA2A)
      -> BruteForceMotifSearch
      ↓
    * Reformulating the "Motif finding problem"
      -> "Equivalent Motif Finding Problem"
         exploring all potential k-mer sonsensus strign frist
         => then find the best possible collection of Motifs for each consensus string
      -> comparing "Motif finding problem" with "Equivalent motif finding problem"
         TODO: What is the key difference?
      -> NEXT: "Median string problem" (BA2B)

Info:
                          T  C  G  G  G  G  g  T  T  T  t  t
                          c  C  G  G  t  G  A  c  T  T  a  C
                          a  C  G  G  G  G  A  T  T  T  t  C
                          T  t  G  G  G  G  A  c  T  T  t  t
    Motifs                a  a  G  G  G  G  A  c  T  T  C  C
                          T  t  G  G  G  G  A  c  T  T  C  C
                          T  C  G  G  G  G  A  T  T  c  a  t
                          T  C  G  G  G  G  A  T  T  c  C  t
                          T  a  G  G  G  G  A  a  c  T  a  C
                          T  C  G  G  G  t  A  T  a  a  C  C
    SCORE(Motifs)         3+ 4+ 0+ 0+ 1+ 1+ 1+ 5+ 2+ 3+ 6+ 4 = 30
    CONSENSUS(Motifs)     T  C  G  G  G  G  A  T  T  T  C  C

Thoughts:
- What exactly is the meaning of this exercise? What are the problems?
  - finding motifs from multiple DNA string (e.g. 15-mer with mutations from each of ten DNA string )
  - Can this exercise be solved by functions in BA1B or BA1I? No.
    cuz ...
       i.   too slow
       ii.  mutations --> BA1B (x)
       iii. ten randomly generated DNA strings
            --> BA1I (x), solution of BA1I is only for single text
            --> concatenate ten DNA strings is indadequate (incorrect model of biological problem)

Terminology:
- motifs (k-mers)
  ; candiates for regulatory motif (or transcription factor binding site)
  ; a collection of k-mers, one from each string in DNA,
    minimizing SCORE(Motifs) among all possible choices of k-mers
- d(Pattern, text)
  ┌──────────────────────────────────────────────────────────────────────┐
  │  d(Pattern, Text) = min          HammingDistance(Pattern, Pattern')  │
  │                  all k-mers Pattern' in Text                         │
  └──────────────────────────────────────────────────────────────────────┘
- d(Pattern, Motifs)
  ┌───────────────────────────────────────────┐
  │                  t                        │
  │  d(Pattern, Dna) = Σ  d(Pattern, Dna_i)   │
  │                  i=0                      │
  └───────────────────────────────────────────┘

Plan 1.
  * brute-force
  * steps:
    a. candidates: all 4^k k-mer patterns as in frequency array (instead of k-mers from a string)
    b. find all (k,d)-mismatched k-mers from the first DNA string
    c. for the rest of DNA strings, do the same thing as in step (b), and intersect two sets

Plan 2.
  * using BA1N neighbors function, set intersection
    a. get all d-neighbors from all k-mers of the first DNA string
    b. repeat set intersection between the first d-neighbors and each of d-neighbor from each DNA string
  * sample dataset
    dnas     k=3    d=1                                      motifs
    -----------------------------------------------------------------
    ATTTGGC  ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT  ATA ATT GTT TTT
             TTT    TTA TTC TTG TAT TCT TGT ATT CTT GTT TTT
             TTG    TTA TTC TAG TCG TGG ATG CTG GTG TTG TTT
             TGG    TGA TGC TAG TCG AGG CGG GGG TGG TTG TGT
             GGC    GGA GAC GCC AGC CGC GGC TGC GTC GGG GGT
    TGCCTTA  TGC    TGA TAC TCC AGC CGC GGC TGC TTC TGG TGT
             GCC    GCA GAC ACC CCC GCC TCC GGC GTC GCG GCT
             CCT    CCA CCC CCG CAT ACT CCT GCT TCT CGT CTT
             CTT    CTA CTC CTG CAT CCT CGT ATT CTT GTT TTT
             TTA    TAA TCA TGA ATA CTA GTA TTA TTC TTG TTT
    CGGTATC  CGG    CGA CGC CAG CCG AGG CGG GGG TGG CTG CGT
             GGT    GGA GGC GGG GAT GCT AGT CGT GGT TGT GTT
             GTA    GAA GCA GGA ATA CTA GTA TTA GTC GTG GTT
             TAT    TAA TAC TAG AAT CAT GAT TAT TCT TGT TTT
             ATC    ATA AAC ACC AGC ATC CTC GTC TTC ATG ATT
    GAAAATT  GAA    AAA CAA GAA TAA GCA GGA GTA GAC GAG GAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAA    AAA CAA GAA TAA ACA AGA ATA AAC AAG AAT
             AAT    AAA AAC AAG AAT CAT GAT TAT ACT AGT ATT
             ATT    ATA ATC ATG AAT ACT AGT ATT CTT GTT TTT

═════════════════════════════════════════════════

References:
- How to check if a map contains a key in Go?
  https://stackoverflow.com/questions/2050391/how-to-check-if-a-map-contains-a-key-in-go
    if val, ok := dict["foo"]; ok {
        //do something here
    }
- How to get intersection of two slice in golang?
  https://stackoverflow.com/questions/44956031/how-to-get-intersection-of-two-slice-in-golang
- Delete key in map
  https://stackoverflow.com/questions/1736014/delete-key-in-map
- Getting a slice of keys from a map
  https://stackoverflow.com/questions/21362950/getting-a-slice-of-keys-from-a-map

************************************************/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"
)

func Neighbors(s string, h int) []string {
	v := []string{s}
	const sym = "A C T G"
	buf := []byte(s)
	var f func([]byte, int)
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

func PatternWithinDistance(Pattern, seq string, d int) bool {
pos:
	for i, j := 0, len(Pattern); j <= len(seq); i, j = i+1, j+1 {
		h := 0
		for i, b := range []byte(seq[i:j]) {
			if Pattern[i] != b {
				if h == d {
					continue pos
				}
				h++
			}
		}
		return true
	}
	return false
}

func ImplantedMotifs(Dna []string, k, d int) []string {
	Patterns := map[string]string{}
	for _, seq := range Dna {
		for i, j := 0, k; j <= len(seq); i, j = i+1, j+1 {
		hv:
			for _, Patternʹ := range Neighbors(seq[i:j], d) {
				for _, seq := range Dna {
					if !PatternWithinDistance(Patternʹ, seq, d) {
						continue hv
					}
				}
				Patterns[string(Patternʹ)] = Patternʹ
			}
		}
	}
	r := make([]string, len(Patterns))
	i := 0
	for _, r[i] = range Patterns {
		i++
	}
	return r
}

func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02a.txt")
	var nums [2]int
	for i, s := range strings.Split(lines[0], " ") {
		nums[i], _ = strconv.Atoi(s)
	}
	k := nums[0]
	d := nums[1]
	dnas := lines[1:]

	// check time
	start := time.Now()
	for _, m := range ImplantedMotifs(dnas, k, d) {
		fmt.Print(m, " ")
	}
	fmt.Println()
	elapsed := time.Since(start)
	fmt.Println("Execution time:", elapsed)

}

/* helper fx: read lines */
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
