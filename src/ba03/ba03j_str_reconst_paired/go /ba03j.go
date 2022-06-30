/*
 */

package main

import (
	"bufio"
	"os"
)

func main() {
	lines := ReadLines("/home/wsl/rosalind/data/ba03a.txt")

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
