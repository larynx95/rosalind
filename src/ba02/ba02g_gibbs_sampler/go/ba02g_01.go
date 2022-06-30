// code by sonia

package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

func RandomKmers(Dna []string, k int) []string {
	kmers := make([]string, len(Dna))
	for i, s := range Dna {
		if len(s) < k {
			return nil
		}
		j := rand.Intn(len(s) - k + 1)
		kmers[i] = s[j : j+k]
	}
	return kmers
}

func Score(kmers []string) (s int) {
	k := len(kmers[0])
	for i := 0; i < k; i++ {
		var f [7]int
		max := 0
		for _, m := range kmers {
			b2 := m[i] & 6
			f[b2]++
			if f[b2] > max {
				max = f[b2]
			}
		}
		s += len(kmers) - max
	}
	return
}

type Profile [][4]float64

func NewLaplaceProfile(kmers []string) Profile {
	inc := 1 / float64(len(kmers))
	p := make(Profile, len(kmers[0]))
	for i := range p {
		p[i] = [4]float64{inc, inc, inc, inc}
	}
	for _, m := range kmers {
		for i, b := range m {
			p[i][b>>1&3^b>>2&1] += inc
		}
	}
	return p
}

func (p Profile) KmerProbability(kmer string) float64 {
	pr := 1.
	for i, b := range []byte(kmer) {
		pr *= p[i][b>>1&3^b>>2&1]
	}
	return pr
}

func RandWeighted(weights []float64) (n int) {
	sum := 0.
	for _, w := range weights {
		sum += w
	}
	f := rand.Float64() * sum
	c := 0.
	for v, w := range weights {
		c += w
		if c >= f {
			return v
		}
	}
	return len(weights) - 1
}

func GibbsSampler(Dna []string, k, N int) (bestMotifs []string, bestScore int) {
	motifs := RandomKmers(Dna, k)
	bestMotifs = make([]string, len(motifs))
	for i, s := range motifs {
		bestMotifs[i] = s
	}
	bestScore = Score(bestMotifs)
	p := make([]float64, len(Dna[0])-k+1)
	for j := 0; j < N; j++ {
		i := rand.Intn(len(Dna))
		motifs[i] = motifs[0]
		pf := NewLaplaceProfile(motifs[1:])
		s := Dna[i]
		p = p[:len(s)-k+1]
		for x := range p {
			p[x] = pf.KmerProbability(s[x : x+k])
		}
		x := RandWeighted(p)
		motifs[i] = s[x : x+k]
		if s := Score(motifs); s < bestScore {
			bestScore = s
			copy(bestMotifs, motifs)
		}
	}
	return
}

func Sample20(Dna []string, k, N int) (motifs []string) {
	motifs, min := GibbsSampler(Dna, k, N)
	for i := 1; i < 20; i++ {
		if m, h := GibbsSampler(Dna, k, N); h < min {
			motifs = m
			min = h
		}
	}
	return
}

func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba02g.txt")
	temp := make([]int, 3)
	for i, str := range strings.Split(lines[0], " ") {
		temp[i], _ = strconv.Atoi(str)
	}
	k := temp[0]
	// t := temp[1]
	n := temp[2]
	dnas := lines[1:]

	for _, m := range Sample20(dnas, k, n) {
		fmt.Println(m)
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
