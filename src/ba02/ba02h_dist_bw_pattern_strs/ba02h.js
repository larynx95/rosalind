/*
Rosalind: BA2H
Implement DistanceBetweenPatternAndStrings

The first potential issue with implementing MedianString from "Find a Median String" is
writing a function to compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai),
the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}.
This task is achieved by the following pseudocode.

  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

Sample Dataset
AAA
TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

Sample Output
5

═════════════════════════════════════════════════

Pseudocode:
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ DISTANCEBETWEENPATTERNANDSTRINGS(Pattern, Dna)                         ║
  ║     k <- |Pattern|                                                     ║
  ║     distance <- 0                                                      ║
  ║     for each string Text in Dna                                        ║
  ║         HammingDistance <- infinite                                    ║
  ║         for each k-mer Pattern’ in Text                                ║
  ║             if HammingDistance > HAMMINGDISTANCE(Pattern, Pattern’)    ║
  ║                 HammingDistance <- HAMMINGDISTANCE(Pattern, Pattern’)  ║
  ║         distance <- distance + HammingDistance                         ║
  ║     return distance                                                    ║
  ╚════════════════════════════════════════════════════════════════════════╝
  ╔════════════════════════════════════════════════════════════════════════╗
  ║ MEDIANSTRING(Dna, k)                                                   ║
  ║     distance 1                                                         ║
  ║     for i 0 to 4^k - 1                                                 ║
  ║         pattern <- NumberToPattern(i, k)                               ║
  ║         if distance > DistanceBetweenPaternAndStrings(Pattern, Dna)    ║
  ║             distance <- DistanceBetweenPatternAndStrings(Patern, Dna)  ║
  ║             Median <- Pattern                                          ║
  ║     return Median                                                      ║
  ╚════════════════════════════════════════════════════════════════════════╝
*/

// BA03A: all k-mers
function* all_kmers_gen(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length - k + 1; i++) yield text.slice(i, i + k);
}

// BA01G: hamming distance
function hamming_distance(astr, bstr) {
    /**
     * (str,str) -> int
     * returns a Hamming distance of two strings (BA1G)
     */
    let distance = 0;
    for (let i = 0; i < astr.length; i++)
        if (astr[i] != bstr[i]) distance++;
    return distance;
}

// BA02H: distance between pattern and strings
function distance_between_pattern_and_strings(pattern, dnas) {
    const k = pattern.length;
    let dist = 0;
    for (const dna of dnas) {
        let min_hdist = Infinity;
        for (const kmer of all_kmers_gen(dna, k)) {
            const hdist = hamming_distance(pattern, kmer);
            if (min_hdist > hdist) min_hdist = hdist;
        }
        dist += min_hdist;
    }
    return dist;
}

// BA01M: number to pattern
function number_to_pattern(num, k) {
    /**
     * (int,int) -> str
     * returns a pattern from a number
     * number_to_pattern(15, 2) == "TT"
     */
    pattern = ""
    for (let i = 0; i < k; i++) {
        temp = Math.floor(num / Math.pow(4, (k - i - 1)));
        if (temp == 0) pattern += 'A';
        else if (temp == 1) pattern += 'C';
        else if (temp == 2) pattern += 'G';
        else if (temp == 3) pattern += 'T';
        num -= temp * Math.pow(4, k - i - 1);
    }
    return pattern
}

// BA02B, BA02H: median string
function median_string(dnas, k) {
    let min_dist = Infinity;
    let medians = new Set();
    for (let i = 0; i < Math.pow(4, k); i++) {
        const pattern = number_to_pattern(i, k);
        const dist = distance_between_pattern_and_strings(pattern, dnas);
        if (dist < min_dist) {
            mind_dist = dist;
            medians = new Set([pattern]);
        } else if (dist == min_dist) {
            medians.add(pattern);
        }
    }
    return medians;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba02h.txt').toString().split("\n");
    const pattern = lines[0];
    const dnas = lines[1].split(/\s/);

    let startTime = performance.now()
    let dist = distance_between_pattern_and_strings(pattern, dnas);
    console.log(dist);
    console.log(`${performance.now() - startTime} milliseconds`)

    startTime = performance.now()
    const mds = median_string(dnas, 3);
    let ans = "";
    for (const m of mds) {
        ans += m + " ";
    }
    console.log(ans);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
