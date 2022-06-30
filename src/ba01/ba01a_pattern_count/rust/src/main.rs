/*
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

Plan 1.
  - using 'math_indices()'
    failed, find 'GCGCG' returns only [0], not [0,2]

Plan 2.
  - get substring, compare, count
    Code is verbose in Rust.
  - chars() + take(), skip() + nth()
  - char_indices() + nth()
  - chars() and enumerate()
  - bytes()

═════════════════════════════════════════════════

References:
- What are the differences between Rust's `String` and `str`?
  https://stackoverflow.com/questions/24158114/what-are-the-differences-between-rusts-string-and-str
- Is there a method like JavaScript's substr in Rust?
  https://stackoverflow.com/questions/37157926/is-there-a-method-like-javascripts-substr-in-rust
- pub fn chars(&self) -> Chars
  String to iterator
  http://web.mit.edu/rust-lang_v1.26.0/arch/amd64_ubuntu1404/share/doc/rust/html/std/primitive.str.html#method.chars
- string length: len(), chars().count()
  'String' in Rust is a UTF-8–encoded, growable string.
  Get the String length in characters in Rust
  https://stackoverflow.com/questions/46290655/get-the-string-length-in-characters-in-rust
- What are the differences between Rust's `String` and `str`?
  https://stackoverflow.com/questions/24158114/what-are-the-differences-between-rusts-string-and-str
  String : mutable, dynamic heap string type,  growable, owned UTF-8 byte array, std::string (C++)
  str    : immutable, string slice, reference to UTF-8 bytes array, fixed length, char* (C++)
- What’s the difference between &String and &str?
  https://users.rust-lang.org/t/whats-the-difference-between-string-and-str/10177
  '&str' can be used for substrings, i.e. they are slices.
  '&String' references always the whole string.
*/

// #!/usr/bin/env rust
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01A: get substring, compare, count
fn count_pattern(text: String, pattern: String) -> u8 {
    let length = text.len(); // len()
    let k = pattern.chars().count(); // converted to iterator, then counted
    let mut count = 0; // mutable variable, default integer type is i32
    for (i, _) in text.char_indices().take(length - k + 1) {
        // a style of iterator in for-loop
        let subs = &text[i..(i + k)]; // :: &str
        if subs.to_string() == pattern {
            count += 1;
        }
    }
    count
}

// BA01A: pattern count, chaining - Is this good?
fn count_pattern_chaining(text: &str, pattern: &str) -> usize {
    let k: usize = pattern.len();
    text.char_indices()
        .map(|(i, _)| i)
        .filter(|&i| i < text.len() - k + 1)
        .filter(|&i| &text[i..(i + k)] == pattern)
        .count() // without semicolon --> return
}

fn main() {
    // practice: String, &String, &str
    let s = "Hello".to_string();
    let a = &s; // &String
    let b = &s[0..2];
    println!("{}", s);
    println!("{}", a);
    println!("{}", &s == a);
    println!("{}", b);

    // faild trial
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01a.txt"); // ::Vec<String>
    let text = &lines[0]; //:: &String
    let pattern = &lines[1]; // :: &String
    let indices_wrong = count_patterns_wrong(text, pattern);
    println!("{:?}", indices_wrong);

    // success
    let start = Instant::now();
    let count = count_pattern(text.to_string(), pattern.to_string());
    println!("{}", count);
    println!("Execution time: {:?}", start.elapsed()); // 12.2 us

    // chaining
    let start = Instant::now();
    let count = count_pattern_chaining(text, pattern);
    println!("{}", count);
    println!("Execution time: {:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}

// creat a custom trait, and impl it to String type... Go back to Python!
trait StringUtils {
    fn substring(&self, start: usize, len: usize) -> Self;
}

impl StringUtils for String {
    fn substring(&self, start: usize, len: usize) -> Self {
        self.chars().skip(start).take(len).collect()
    }
}

// using match_indices(), but failed
fn count_patterns_wrong(text: &str, pattern: &str) -> Vec<usize> {
    let indices: Vec<_> = text.match_indices(pattern).map(|(i, _)| i).collect();
    return indices;
}
