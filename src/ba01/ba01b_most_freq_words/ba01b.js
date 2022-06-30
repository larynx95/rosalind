/*
Rosalind: BA1B
Find the Most Frequent Words in a String

We say that Pattern is a most frequent k-mer in Text
if it maximizes Count(Text, Pattern) among all k-mers.
For example, "ACTAT" is a most frequent 5-mer
in "ACAACTATGCATCACTATCGGGAACTATCCT",
and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4

Sample Output
CATG GCAT

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> HERE: find frequent words in a string (BA1B)
      -> NEXT: find all occurrence of a pattern in a string (BA1D)

Plan 1.
- writing down
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
- Pseudocode:
╔═════════════════════════════════════════════════╗
║ COMPUTINGFREQUENCIES(Text, k)                   ║
║   for i <- 0 to 4^k- 1                          ║
║     FREQUENCYARRAY(i) <- 0                      ║
║   for i <- 0 to |Text| - k                      ║
║     Pattern <- Text(i, k)                       ║
║     j <- PATTERNTONUMBER(Pattern)               ║
║     FREQUENCYARRAY(j) <- FREQUENCYARRAY(j) + 1  ║
║   return FREQUENCYARRAY                         ║
╚═════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════╗
║ FREQUENTWORDS(Text, k)                          ║
║   FrequentPatterns an empty set                 ║
║     for i <- 0 to |Text| - k                    ║
║     Pattern <- the k-mer Text(i, k)             ║
║     COUNT(i) <- PATTERNCOUNT(Text, Pattern)     ║
║   maxCount maximum value in array COUNT         ║
║   for i <- 0 to |Text| - k                      ║
║     if COUNT(i) = maxCount                      ║
║       add Text(i, k) to FrequentPatterns        ║
║   remove duplicates from FrequentPatterns       ║
║   return FrequentPatterns                       ║
╚═════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════╗
║ FASTERFREQUENTWORDS(Text , k)                     ║
║   FrequentPatterns <- an empty set                ║
║   FREQUENCYARRAY <- COMPUTINGFREQUENCIES(Text, k) ║
║   maxCount <- maximal value in FREQUENCYARRAY     ║
║   for i <- 0 to 4k - 1                            ║
║     if FREQUENCYARRAY(i) = maxCount               ║
║       Pattern <- NUMBERTOPATTERN(i, k)            ║
║       add Pattern to the set FrequentPatterns     ║
║   return FrequentPatterns                         ║
╚═══════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════╗
║ FINDINGFREQUENTWORDSBYSORTING(Text , k)              ║
║   FrequentPatterns <- an empty set                   ║
║   for i <- 0 to |Text| - k                           ║
║     Pattern <- Text(i, k)                            ║
║     INDEX(i) <- PATTERNTONUMBER(Pattern)             ║
║     COUNT(i) <- 1                                    ║
║   SORTEDINDEX <- SORT(INDEX)                         ║
║   for i <- 1 to |Text| - k                           ║
║     if SORTEDINDEX(i) = SORTEDINDEX(i - 1)           ║
║       COUNT(i) = COUNT(i - 1) + 1                    ║
║   maxCount <- maximum value in the array COUNT       ║
║   for i <- 0 to |Text| - k                           ║
║     if COUNT(i) = maxCount                           ║
║       Pattern <- NUMBERTOPATTERN(SORTEDINDEX(i), k)  ║
║       add Pattern to the set FrequentPatterns        ║
║   return FrequentPatterns                            ║
╚══════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- How to count string occurrence in string?
  https://stackoverflow.com/questions/4009756/how-to-count-string-occurrence-in-string
- Checking if a key exists in a JavaScript object?
  https://stackoverflow.com/questions/1098040/checking-if-a-key-exists-in-a-javascript-object
- javascript set
  https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Set
- Iterate over set elements
  https://stackoverflow.com/questions/16401216/iterate-over-set-elements
- How to perform an integer division, and separately get the remainder, in JavaScript?
  https://stackoverflow.com/questions/4228356/how-to-perform-an-integer-division-and-separately-get-the-remainder-in-javascr
- Most efficient way to create a zero filled JavaScript array?
  https://stackoverflow.com/questions/1295584/most-efficient-way-to-create-a-zero-filled-javascript-array
- Is there a “not in” operator in JavaScript for checking object properties?
  https://stackoverflow.com/questions/7972446/is-there-a-not-in-operator-in-javascript-for-checking-object-properties
- Find the min/max element of an array in JavaScript
  https://stackoverflow.com/questions/1669190/find-the-min-max-element-of-an-array-in-javascript
*/

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

// BA01M: number to pattern
function number_to_pattern_rec(num, k) {
    /**
     * (int,int) -> str
     * returns a pattern from a number
     * number_to_pattern_rec(15, 2) == "TT"
     */
    const dic = { 0: "A", 1: "C", 2: "G", 3: "T" };
    if (k == 0) return ""
    else {
        temp = Math.floor(num / Math.pow(4, (k - 1)));
        num -= temp * Math.pow(4, k - 1);
        return dic[temp] + number_to_pattern(num, k - 1)
    }
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

// BA01L: pattern to number
function pattern_to_number_rec(pattern) {
    /**
     * str -> int
     * returns a number from a pattern
     * pattern_to_number_rec("TT") == 15
     */
    const dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
    if (pattern.length == 0) {
        return 0;
    } else {
        return dic[pattern[0]] * Math.pow(4, pattern.length - 1) + pattern_to_number(pattern.slice(1));
    }
}

// BA01K: compute frequency
function compute_freq(text, k) {
    /**
     * (str,int) -> [int]
     * returns a frequency array
     * compute_freq('ACGCGGCTCTGAAA', 2)
     * == [2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]
     */
    let freq_arr = Array(Math.pow(4, k));  // TODO: find the fastest way
    for (let i = 0; i < Math.pow(4, k); i++) freq_arr[i] = 0;
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        freq_arr[num] += 1;
    }
    return freq_arr;
}

// BA01K: compute frequency
function compute_freq_dic(text, k) {
    /**
     * (str,int) -> ({int:int},int)
     * returns a frequency array
     */
    let freq_dic = {};
    let max_cnt = 0;
    for (let i = 0; i < text.length - k + 1; i++) {
        let kmer = text.slice(i, i + k);
        let num = pattern_to_number(kmer);
        if (!(num in freq_dic)) freq_dic[num] = 1;
        else {
            freq_dic[num] += 1;
            if (freq_dic[num] > max_cnt) max_cnt = freq_dic[num];
        }
    }
    return [freq_dic, max_cnt];
}

// BA01B: most frequent words
function most_freq_words(text, k) {
    /**
     * (str,int) -> {str}
     * returns a set of most frequent words
     * dictionary (object version)
     */
    let dict = {}
    let max_cnt = 0
    for (let i = 0; i < text.length - k + 1; i++) {
        subs = text.slice(i, i + k)
        if (subs in dict) {
            let cnt = ++dict[subs];
            if (cnt > max_cnt) max_cnt = cnt;
        } else {
            dict[subs] = 1;
        }
    }
    answer = new Set();
    for (key in dict) {
        if (dict[key] == max_cnt) {
            answer.add(key);
        }
    }
    return answer;
}

// BA01B: most frequent words
function most_freq_words_faster(text, k) {
    /**
     * (str,int) -> {str}
     * returns a set of most frequent words
     */
    let freq_arr = compute_freq(text, k);
    let answer = new Set();
    let max_cnt = Math.max(...freq_arr);
    for (let i = 0; i < freq_arr.length; i++) {
        if (freq_arr[i] == max_cnt) answer.add(number_to_pattern(i, k));
    }
    return answer;
}

// BA01B: most frequent words
function most_freq_words_faster_dic(text, k) {
    /**
     * (str,int) -> {str}
     * returns a set of most frequent words
     */
    let [freq_dic, max_cnt] = compute_freq_dic(text, k);
    let answer = new Set();
    for (const key in freq_dic) {
        if (freq_dic[key] == max_cnt) answer.add(number_to_pattern(key, k));
    }
    return answer;
}

// main
function main() {
    const fs = require('fs');
    let lines = fs.readFileSync('/home/wsl/rosalind/data/ba01b.txt').toString().split("\n");
    const text = lines[0];
    const k = parseInt(lines[1]);

    let startTime = performance.now();
    let words = most_freq_words_faster_dic(text, k);
    let str_answer = "";
    for (let answer of words) {
        str_answer += answer + " ";
    }
    console.log(str_answer);
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
