/*
*/

use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03i.txt");

    // check time
    let start = Instant::now();
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

