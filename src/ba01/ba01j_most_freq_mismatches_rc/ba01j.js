/*
Rosalind: BA1J
Find Frequent Words with Mismatches and Reverse Complements

We now extend “Find the Most Frequent Words with Mismatches in a String”
to find frequent words with both mismatches and reverse complements.
Recall that rPattern refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
Find the most frequent k-mers (with mismatches and reverse complements)
in a DNA string.

Given: A DNA string Text as well as integers k and d.

Return:
All k-mers Pattern
maximizing the sum Count_d(Text, Pattern) + Count_d(Text, rPattern)
over all possible k-mers.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
ATGT ACAT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> find frequent words in a string (BA1B)
      -> find all occurrence of a pattern in a string (BA1D)
      -> Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> hamming distance (BA1G)
      -> Find All Approximate Occurrences of a Pattern in a String (BA1H)
      -> Find the Most Frequent Words with Mismatches in a String (BA1I)
      -> HERE: Find Frequent Words with Mismatches and Reverse Complements (BA1J)
        -> Find the Reverse Complement of a String (BA1C)
      ↓
    * Solving the most frequent words with mismatch by improving algorithms
      -> Generate the Frequency Array of a String (BA1K)
        -> NEXT: pattern to number (BA1L), number to pattern (BA1M)

*/

// BA01C: reverse complement
function reverse_complement(pattern) {
    /**
     * str -> str
     * returns a reverse complement
     * reverse_complement("ACGGTT") == AACCGT
     */
    dic = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' };
    pattern = [...pattern].reverse().join("");
    answer = "";
    for (const chr of pattern) {
        answer += dic[chr];
    }
    return answer;
}

// BA01G: hamming distance
function hamming_distance(astr, bstr) {
    /**
     * (str,str) -> int
     * returns a Hamming distance of two strings
     */
    let distance = 0;
    for (let i = 0; i < astr.length; i++)
        if (astr[i] != bstr[i]) distance++;
    return distance;
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

// BA01L: pattern to number
function pattern_to_number(pattern) {
    /**
     * str -> int
     * returns a number from a pattern
     * pattern_to_number("TT") == 15
     */
    const dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
    let number = 0;
    for (let i = 0; i < pattern.length; i++) {
        let pw = pattern.length - i - 1;
        let nuc = 0;
        let chr = pattern[i];
        if (chr == 'A') nuc = 0;
        else if (chr == 'C') nuc = 1;
        else if (chr == 'G') nuc = 2;
        else if (chr == 'T') nuc = 3;
        let pw_val = Math.pow(4, pw);
        number += pw_val * nuc;
    }
    return number;
}

// BA01N: neighbors
function neighbors_iter(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    let neighborhood = new Set([pattern]);
    for (let i = 0; i < d; i++) {
        for (const pat of neighborhood) {
            neighborhood = new Set([...immediate_neighbors(pat), ...neighborhood]);
        }
    }
    return neighborhood;
}

// BA01N: neighbors
function neighbors(pattern, d) {
    /**
     * (str,int) -> {str}
     * returns a set of d-neighbors
     */
    //
    if (d == 0) return new Set([pattern]);
    if (pattern.length == 0) return new Set();
    if (pattern.length == 1) return new Set(['A', 'C', 'G', 'T']);
    let neighborhood = new Set();
    let suffixneighbors = neighbors(pattern.slice(1), d);
    for (const text of suffixneighbors) {
        if (hamming_distance(text, pattern.slice(1)) < d) {
            for (const nuc of ['A', 'C', 'G', 'T']) {
                neighborhood.add(nuc + text);
            }
        } else {
            neighborhood.add(pattern[0] + text);
        }
    }
    return neighborhood;
}

// BA01I: frequent words mismatches
function freq_words_mismatches(text, k, d) {
    /**
     * (str,int,int) -> {str}
     * returns a set of frequent words with mismatches
     */
    let freq_patterns = new Set();
    let count = {};
    let max_cnt = -1;
    for (let i = 0; i < text.length - k + 1; i++) {
        let neighborhood = neighbors(text.slice(i, i + k), d);
        for (const word of neighborhood) {
            let idx = pattern_to_number(word);
            if (!(idx in count)) count[idx] = 1;
            else count[idx]++;
            if (count[idx] > max_cnt) max_cnt = count[idx];
        }
    }
    for (const num in count) {
        if (count[num] == max_cnt) {
            freq_patterns.add(number_to_pattern(num, k));
        }
    }
    return freq_patterns;
}

// BA01J: frequent words mismatches reverse complement
function freq_words_mismatches_rev_complement(text, k, d) {
    /**
     * (str,int,int) -> {str}
     */
    // create a frequency dictionary
    let count = {};
    for (let i = 0; i < text.length - k + 1; i++) {
        let neighborhood = neighbors(text.slice(i, i + k), d);
        for (const word of neighborhood) {
            let idx = pattern_to_number(word);
            if (!(idx in count)) count[idx] = 1;
            else count[idx]++;
        }
    }
    // find frequent words with mismatches and reverse complement
    // searching and counting simultaneously
    let max_cnt = -1;
    let freq_patterns = new Set();
    for (const num in count) {
        let num_rc = pattern_to_number(reverse_complement(number_to_pattern(num, k)));
        let cnt = count[num] + ((num_rc in count) ? count[num_rc] : 0);  // trim undefined cases
        if (cnt > max_cnt) {
            freq_patterns = new Set([number_to_pattern(num, k)]);
            max_cnt = cnt;
        } else if (cnt == max_cnt) {
            freq_patterns.add(number_to_pattern(num, k));
        }
    }
    return freq_patterns;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba01j.txt').toString().split("\n");
    const text = lines[0];
    const [k, d] = lines[1].split(/\s/).map(Number);

    // check time
    const startTime = performance.now();
    let answer = "";
    for (const word of freq_words_mismatches_rev_complement(text, k, d)) {
        answer += word + " ";
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
