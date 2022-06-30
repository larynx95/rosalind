/*
Rosalind: BA1A
Compute the Number of Times a Pattern Appears in a Text

Rosalind: BA1A
Compute the Number of Times a Pattern Appears in a Text

This is the first problem in a collection of "code challenges"
to accompany Bioinformatics Algorithms:
An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

A k-mer is a string of length k.
We define Count(Text, Pattern) as the number of times
that a k-mer Pattern appears as a substring of Text.
For example,

Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.

We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2)
since we should account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern),
our plan is to "slide a window" down Text,
checking whether each k-mer substring of Text matches Pattern.
We will therefore refer to the k-mer starting at position i of Text as Text(i,k).
Throughout this book, we will often use 0-based indexing,
meaning that we count starting at 0 instead of 1.
In this case, Text begins at position 0 and ends at position |Text| - 1
(|Text| denotes the number of symbols in Text).
For example, if Text = GACCATACTG, then Text(4, 3) = ATA.
Note that the last k-mer of Text begins at position |Text| - k, e.g.,
the last 3-mer of GACCATACTG starts at position 10 - 3 = 7.
This discussion results in the following pseudocode for computing Count(Text, Pattern).

╔═════════════════════════════════════════╗
║ PatternCount(Text, Pattern)             ║
║     count <- 0                          ║
║     for i <- 0 to |Text| - |Pattern|    ║
║         if Text(i, |Pattern|) = Pattern ║
║             count <- count + 1          ║
║     return count                        ║
╚═════════════════════════════════════════╝

Implement PatternCount
Given: {DNA strings}} Text and Pattern.
Return: Count(Text, Pattern).

Sample Dataset
GCGCG
GCG

Sample Output
2

═════════════════════════════════════════════════

References:
- Is there a JavaScript strcmp()?
  https://stackoverflow.com/questions/1179366/is-there-a-javascript-strcmp
- node.js: read a text file into an array. (Each line an item in the array.)
  https://stackoverflow.com/questions/6831918/node-js-read-a-text-file-into-an-array-each-line-an-item-in-the-array
- In Node.js, how do I "include" functions from my other files?
  https://stackoverflow.com/questions/5797852/in-node-js-how-do-i-include-functions-from-my-other-files
*/

// BA01: pattern count
function patternCount_substring1(text, pattern) {
    let count = 0
    for (let i = 0; i < text.length - pattern.length + 1; i++) {
        substring = text.substring(i, i + pattern.length);
        if (substring === pattern) {
            count++;
        }
    }
    return count
}

// BA01: pattern count
function patternCount_substring2(text, pattern) {
    let count = 0
    for (let i = 0; i < text.length - pattern.length + 1; i++) {
        if (text.substring(i, pattern.length).localeCompare(pattern) == 1) {
            count++;
        }
    }
    return count
}

// BA01: pattern count
function patternCount_slice(text, pattern) {
    let count = 0
    for (let i = 0; i < text.length - pattern.length + 1; i++) {
        substring = text.slice(i, i + pattern.length);
        //console.log(typeof substring);
        //console.log(typeof pattern);
        //console.log(substring == pattern, substring, pattern)
        if (substring == pattern) {
            count++;
        }
    }
    return count
}

// main
function main() {
    const fs = require('fs');
    let lines = fs.readFileSync('/home/wsl/rosalind/data/ba01a.txt').toString().split("\n");
    text = lines[0]
    pattern = lines[1]

    // check time
    let startTime = performance.now()
    answer = patternCount_slice(text, pattern)
    console.log(answer)
    let endTime = performance.now()
    console.log(`${endTime - startTime} milliseconds`)
}

// execute main function
main()
