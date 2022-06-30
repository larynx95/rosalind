/*
Rosalind: BA3E (difficulty: 1/5)
Construct the "De Bruijn Graph" of a Collection of k-mers

Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times),
we define "CompositionGraph(Patterns)" as a graph with |Patterns| isolated edges.
Every edge is labeled by a k-mer from Patterns,
and the starting and ending nodes of an edge are labeled
by the prefix and suffix of the k-mer labeling that edge.
We then define the "de Bruijn graph" of Patterns, denoted "DeBruijn(Patterns)",
by gluing identically labeled nodes in CompositionGraph(Patterns),
which yields the following algorithm.

    DEBRUIJN(Patterns)
        represent every k-mer in Patterns as an isolated edge between its prefix and suffix
        glue all nodes with identical labels, yielding the graph DeBruijn(Patterns)
        return DeBruijn(Patterns)

"De Bruijn Graph" from k-mers Problem
Construct the "de Bruijn graph" from a collection of k-mers.

Given: A collection of k-mers Patterns.

Return: The "de Bruijn graph" "DeBruijn(Patterns)", in the form of an adjacency list.

Sample Dataset
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG

Sample Output
AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> PREV: Construct De Bruijn Graph (BA3D)
      -> HERE: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> NEXT: Eulerian Cycle (BA3F)

═════════════════════════════════════════════════

References:
-
*/

use std::collections::HashMap;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03E: De Bruijn graph from k-mer patterns
fn debruijn_from_kmers(patterns: &Vec<String>) -> HashMap<String, Vec<String>> {
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    for pattern in patterns {
        let k = &pattern.chars().count();
        let prefix = &pattern[..(k - 1)];
        let suffix = &pattern[1..];
        if map.contains_key(prefix) {
            map.get_mut(prefix).unwrap().push(suffix.to_string());
        } else {
            map.insert(prefix.to_string(), vec![suffix.to_string()]);
        }
    }
    return map;
}

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03e.txt");

    // check time
    let start = Instant::now();
    for (kmer, values) in debruijn_from_kmers(&lines) {
        print!("{kmer} -> ");
        print!("{}", values.into_iter().collect::<Vec<String>>().join(","));
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
