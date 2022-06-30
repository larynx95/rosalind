/*
Rosalind: BA3C
Construct the Overlap Graph of a Collection of k-mers

In this chapter, we use the terms prefix and suffix
to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.

Given an arbitrary collection of k-mers Patterns,
we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern'
by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern').
The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

Overlap Graph Problem
Construct the overlap graph of a collection of k-mers.

Given:
A collection Patterns of k-mers.

Return:
The overlap graph Overlap(Patterns), in the form of an adjacency list.

Sample Dataset
ATGCG
GCATG
CATGC
AGGCA
GGCAT

Sample Output (1:1)
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * HERE: Construct OverapGraph (BA3C)
      -> NEXT: Reconstruct a String from its k-mer Composition (BA3H)

Plan 1.
  * Questions:
    - functino arg: pass-by-value vs. pass-by-reference
    - iteration of vector vs. &vector

═════════════════════════════════════════════════

References:
* iteration of vector vs. &vector
  - Why can I iterate over a slice twice, but not a vector?
    https://stackoverflow.com/questions/34572784/why-can-i-iterate-over-a-slice-twice-but-not-a-vector
* &String -> String conversion
  - `to_owned()` method
    https://stackoverflow.com/questions/45275362/how-do-i-copy-a-string-from-a-string-while-iterating-through-a-vector
*/

use std::collections::HashMap;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03C: overlap graph
fn overlap_graph(patterns: &Vec<String>) -> HashMap<String, Vec<String>> {
    let mut dic: HashMap<String, Vec<String>> = HashMap::new(); //
    for source in patterns {
        // 'patterns' is '&Vec<String'(reference), not 'Vec<String>'(value)
        for target in patterns {
            if &source[1..] == &target[..(target.len() - 1)] {
                if dic.contains_key(&source.to_string()) {
                    dic.get_mut(source) // get mutable value by key
                        .unwrap() // 'get_mut()` returns Option
                        .push(target.as_str().to_string()); // 'target' is '&String', -> &str -> String
                                                            // or use ;to_owned()` method
                } else {
                    dic.insert(source.to_string(), vec![target.to_string()]);
                }
            }
        }
    }
    return dic;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03c.txt");

    // check time
    let start = Instant::now();
    for (source, targets) in overlap_graph(&lines) {
        print!("{} -> ", source);
        for target in targets {
            print!("{} ", target);
        }
        println!()
    }
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
