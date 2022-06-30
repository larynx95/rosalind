/*
Rosalind: BA1F
Find a Position in a Genome Minimizing the Skew

Define the skew of a DNA string Genome, denoted Skew(Genome),
as the difference between the total number of occurrences of 'G' and 'C' in Genome.
Let Prefixi(Genome) denote the prefix (i.e., initial substring) of Genome of length i.
For example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

C  A  T  G G G C A T C G G C C A T A C  G C  C
0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i
(from 0 to |Genome|).

Sample Dataset
CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG

Sample Output
53 97

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> HERE: minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * NEXT: Most frequent words problem

Plan 1.
  * sample dataset
      C A T G G G C A T C G G C C A T A C G C C
     -1    +1+1+1-1    -1+1+1-1-1      -1+1-1-1
    0-1-1-1 0 1 2 1 1 1 0 1 2 1 0 0 0 0-1 0-1-2
*/

use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA01F: min skewness
fn min_skewness(text: &str) -> Vec<i32> {
    let mut skewness = 0; // skewness, mutable
    let mut indices: Vec<i32> = vec![0]; // result, mutable
    let mut min_skew = i32::MAX;  // minimum skewness
    for (i, ch) in text.char_indices() {
        skewness += match ch {
            'C' => -1,  // need i32, not usize
            'G' => 1,
            _ => 0,
        };
        if skewness < min_skew {
            min_skew = skewness;
            indices = vec![i as i32 + 1];
        } else if skewness == min_skew {
            indices.push(i as i32 + 1);
        }
    }
    return indices;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba01f.txt");
    let text = &lines[0];

    let start = Instant::now();
    let answer: Vec<i32> = min_skewness(text);
    for val in answer {
        print!("{} ", val);
    }
    println!();
    println!("Execution time: {:?}", start.elapsed());
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    // file automatically dropped when it is out-of scope.
    // no need to 'close' file
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
