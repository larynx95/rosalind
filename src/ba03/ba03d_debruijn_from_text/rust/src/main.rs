/*
Rosalind: BA3D
Construct the "De Bruijn[broin]" Graph of a String

Given a genome Text, PathGraph_k(Text) is the path consisting of |Text| - k + 1 edges,
where the i-th edge of this path is labeled by the i-th k-mer in Text
and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text.
The de Bruijn graph DeBruijn_k(Text) is formed
by gluing identically labeled nodes in PathGraph_k(Text).

De Bruijn Graph from a String Problem
Construct the de Bruijn graph of a string.

Given: An integer k and a string Text.

Return:DeBruijn_k(Text), in the form of an adjacency list.

Sample Dataset
4
AAGATTCTCTAC

Sample Output (1:many)
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC       <-- TODO: Is duplication allowed?
TTC -> TCT

═════════════════════════════════════════════════

Info.
  * what if there're duplication?
  (1)   (2)
  AAGATTAAGATT, 4
  AAGA:  AAG -> AGA
  AGAT:  AGA -> GAT
  GATT:  GAT -> ATT
  ATTA:  ATT -> TTA
  TTAA:  TTA -> TAA
  TAAG:  TAA -> AAG
  AAGA:  AAG -> AGA
  AGAT:  AGA -> GAT
  GATT:  GAT -> ATT

═════════════════════════════════════════════════

References:
- How to join elements of HashSet into a String with a delimiter
  https://stackoverflow.com/questions/47578011/how-to-join-elements-of-hashset-into-a-string-with-a-delimiter
*/

use std::collections::HashMap;
use std::collections::HashSet;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// BA03D: De Bruijn graph from a text, allowing duplicates
fn debruijn_from_text(text: &str, k: usize) -> HashMap<String, Vec<String>> {
    // used HashSet to avoid duplicates
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..(text.chars().count() - k + 1) {
        let prefix = &text[i..(i + k - 1)];
        let suffix = &text[(i + 1)..(i + k)];
        if !map.contains_key(prefix) {
            // if key not in map
            map.insert(prefix.to_string(), vec![suffix.to_string()]);
        } else {
            map.get_mut(prefix).unwrap().push(suffix.to_string());
        }
    }
    return map;
}

// BA03D: De Bruijn graph from a text, no duplicate version
fn debruijn_from_text_hashset(text: &str, k: usize) -> HashMap<String, HashSet<String>> {
    // used HashSet to avoid duplicates
    let mut map: HashMap<String, HashSet<String>> = HashMap::new();
    for i in 0..(text.chars().count() - k + 1) {
        let prefix = &text[i..(i + k - 1)];
        let suffix = &text[(i + 1)..(i + k)];
        if !map.contains_key(prefix) {
            map.insert(prefix.to_string(), HashSet::from([suffix.to_string()]));
        } else {
            map.get_mut(prefix).unwrap().insert(suffix.to_string());
        }
    }
    return map;
}
/*
TODO: Which is correct?
text: "AAGATTCTCTACAAGATTCTCTAC" ("AAGATTCTCTAC" + "AAGATTCTCTAC")
k: 4
the first function                the second function
-----------------------------------------------------
AAG -> AGA,AGA                    TTC -> TCT
ATT -> TTC,TTC                    CTA -> TAC
TTC -> TCT,TCT                    TAC -> ACA
TAC -> ACA                        ATT -> TTC
CAA -> AAG                        TCT -> CTC,CTA
AGA -> GAT,GAT                    CTC -> TCT
GAT -> ATT,ATT                    ACA -> CAA
TCT -> CTC,CTA,CTC,CTA            CAA -> AAG
CTC -> TCT,TCT                    AAG -> AGA
CTA -> TAC,TAC                    AGA -> GAT
ACA -> CAA                        GAT -> ATT
*/

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03d.txt");
    let k: usize = lines[0].parse::<usize>().unwrap();
    let text = &lines[1];

    // check time
    let start = Instant::now();
    for (key, values) in debruijn_from_text(text, k) {
        print!("{key} -> ");
        print!("{}", values.into_iter().collect::<Vec<String>>().join(","));
        println!()
    }
    println!("{:?}\n", start.elapsed());

    // check duplication
    println!("result with duplication==============");
    for (key, values) in debruijn_from_text("AAGATTCTCTACAAGATTCTCTAC", k) {
        print!("{key} -> ");
        print!("{}", values.into_iter().collect::<Vec<String>>().join(","));
        println!()
    }

    println!("result without duplication===========");
    for (key, values) in debruijn_from_text_hashset("AAGATTCTCTACAAGATTCTCTAC", k) {
        print!("{key} -> ");
        print!("{}", values.into_iter().collect::<Vec<String>>().join(","));
        println!()
    }
}

// helper function: read a text file, create a Vector of strings
fn lines_from_file(filename: &str) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
