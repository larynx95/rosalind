/*
Rosalind: BA3B
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern[1], ... , Pattern[n]
such that the last k - 1 symbols of Pattern[i] are equal
to the first k - 1 symbols of Pattern[i+1] for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> PREV: Generate the k-mer Composition of a String (BA3A)
      -> HERE: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * NEXT: Construct OverapGraph (BA3C)

Plan 1.
  * start from the first pattern, add the last characters from each pattern
    ACCGA
     CCGAA
      CGAAG
       GAAGC
        AAGCT

═════════════════════════════════════════════════
References:
- Get the String length in characters in Rust
  https://stackoverflow.com/questions/46290655/get-the-string-length-in-characters-in-rust
- How do I concatenate strings?
  https://stackoverflow.com/questions/30154541/how-do-i-concatenate-strings
*/

use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03B: string spelled by genome path
fn str_spelled_by_genome_path(patterns: &Vec<String>) -> String {
    let mut result = patterns[0].clone();  // TODO: find the better way.
    for pattern in &patterns[1..] {
        result.push(pattern.chars().last().unwrap()); // last character
    }
    return result;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03b.txt");

    // check time
    let start = Instant::now();
    println!("{}", str_spelled_by_genome_path(&lines));
    println!("{:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
