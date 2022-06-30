/*
Rosalind: BA3H
Reconstruct a String from its k-mer Composition

String Reconstruction Problem
Reconstruct a string from its k-mer composition.

Given: An integer k followed by a list of k-mers Patterns.

Return: A string Text with k-mer composition equal to Patterns.
(If multiple answers exist, you may return any one.)

Sample Dataset
4
CTTA
ACCA
TACC
GGCT
GCTT
TTAC

Sample Output
GGCTTACCA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * PREV: Construct OverapGraph (BA3C)
      -> HERE: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> NEXT: Construct De Bruijn Graph (BA3D)

Plan 1.
  * multiple strings -> overlap graph -> Hamiltonian path
  * split k-mer to two (k-1)-mers (prefix and suffix)
    CTTA : CTT,TTA
    ACCA : ACC,CCA
    TACC : TAC,ACC
    GGCT : GGC,GCT
    GCTT : GCT,CTT
    TTAC : TTA,TAC
    prefixes: [CCT,ACC,TAC,GGC,GCT,TTA]
    suffixes: [TTA,CCA,ACC,GCT,CTT,TAC]

Plan 2.
  * using Eulerian path (BA3G) algorithm
  * multiple strings -> De Bruijn graph -> Eulerian cycle -> Eulerian path

═════════════════════════════════════════════════

References:
-
*/

use rand::Rng;
use std::collections::HashMap;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03h.txt");

    // check time
    let start = Instant::now();
    let mut graph = debruijn_from_kmers(&lines[1..lines.len()].to_vec());
    // print_graph(&graph);

    let path = eulerian_path(&mut graph);
    let mut answer = path[0].clone();
    for node in path[1..path.len()].iter() {
        answer.push(node.chars().last().unwrap());
    }
    println!("{}", answer);
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

// helper fx: read graph
fn read_graph(lines: &Vec<String>) -> HashMap<String, Vec<String>> {
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();
    for line in lines {
        let splitted = line.split(" -> ").collect::<Vec<&str>>();
        let key: String = splitted[0].to_string();
        let mut values: Vec<String> = vec![];
        for value in splitted[1].split(',') {
            values.push(value.to_string());
        }
        graph.insert(key, values);
    }
    return graph;
}

// helper fx: print graph
fn print_graph(graph: &HashMap<String, Vec<String>>) {
    for (key, values) in graph {
        print!("{key} -> ");
        println!("{}", values.join(","));
    }
}

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

// walk
fn walk(graph: &mut HashMap<String, Vec<String>>, source: &mut String) -> Vec<String> {
    let mut cycle: Vec<String> = vec![source.to_string()];
    // if target nod exists, add target node to cycle repeatedly
    while let Some(nodes) = graph.get(source) {
        // if value (Vec<String) is empty, break while loop
        if nodes.is_empty() {
            break;
        }
        // get target node, and remove it from graph[source] - Rust's pop() function is very useful.
        let target = graph.get_mut(source).unwrap().pop().unwrap();
        // add target to cycle
        cycle.push(target.clone());
        // update source - The type of source is &String, so life-time doesn't matter.
        *source = target;
    }
    return cycle;
}

// is_depleted
fn is_depleted(graph: &mut HashMap<String, Vec<String>>) -> bool {
    for values in graph.values_mut() {
        if !values.is_empty() {
            return false;
        }
    }
    return true;
}

// BA03F: Eulerian cycle, destructive
fn eulerian_cycle(graph: &mut HashMap<String, Vec<String>>) -> Vec<String> {
    // select start node randomly
    let mut keys: Vec<String> = vec![];
    for (key, _) in graph.iter() {
        keys.push(key.to_string());
    }
    let idx = rand::thread_rng().gen_range(0..keys.len());
    let start = &mut keys[idx];
    // get the first cycle (will be updated repeatedly)
    let mut cycle = walk(graph, start);
    // loop if not all the values are empty vector
    while !is_depleted(graph) {
        let mut new_cycle: Vec<String> = vec![];
        // select another start node, get new cycle (Don't make temporary variable. - Value does not live long enough)
        for node in cycle.iter_mut() {
            /*
            // Be careful! 'graph.get(node)' can be None
            if !graph.get(node).unwrap().is_empty() {
                new_cycle = walk(graph, node);
                break;
            }
            */
            /*
            // Be careful! 'vec![]' is Some(vec![]), not None
            if let Some(_) = graph.get(node) {
                new_cycle = walk(graph, node);
                break;
            }
            */
            if let Some(v) = graph.get(node) {
                if !v.is_empty() {
                    new_cycle = walk(graph, node);
                    break;
                }
            }
        }
        // get index
        let mut idx: usize = 0;
        for (i, elem) in cycle.iter().enumerate() {
            if elem == &new_cycle[0] {
                idx = i;
            }
        }
        // merge new_cycle to main cycle
        let fst_part: Vec<String> = cycle[..idx].to_vec();
        let lst_part: Vec<String> = cycle[idx + 1..].to_vec();
        cycle = fst_part;
        cycle.extend(new_cycle);
        cycle.extend(lst_part)
    }
    return cycle;
}

// get degrees
fn get_degrees(graph: &HashMap<String, Vec<String>>) -> HashMap<String, (usize, usize)> {
    let mut degrees: HashMap<String, (usize, usize)> = HashMap::new();
    for (key, values) in graph.iter() {
        if degrees.contains_key(key) {
            degrees.get_mut(key).unwrap().1 += values.len();
        } else {
            degrees.insert(key.to_string(), (0, values.len()));
        }
        for val in values {
            if degrees.contains_key(val) {
                degrees.get_mut(val).unwrap().0 += 1;
            } else {
                degrees.insert(val.to_string(), (1, 0));
            }
        }
    }
    return degrees;
}

// find start, sink node
fn find_start_sink(degrees: &HashMap<String, (usize, usize)>) -> (String, String) {
    let mut start: String = String::new();
    let mut sink: String = String::new();
    for (key, tup) in degrees {
        if tup.0 < tup.1 {
            start = key.to_string();
        } else if tup.0 > tup.1 {
            sink = key.to_string();
        }
    }
    return (start, sink);
}

// BA03G: Eulerian path
fn eulerian_path(graph: &mut HashMap<String, Vec<String>>) -> Vec<String> {
    // get degrees, find start, sink node
    let degrees = get_degrees(graph);
    let tup = find_start_sink(&degrees);
    let sink = &tup.1;

    // connect sink to start
    graph.insert(tup.1.clone(), vec![tup.0]);

    // get cycle
    let cycle = eulerian_cycle(graph);

    // find index of sink node in cycle
    let idx = cycle
        .iter()
        .enumerate()
        .filter(|&(_, e)| e == sink)
        .collect::<Vec<(usize, &String)>>()
        .pop()
        .unwrap()
        .0;

    // if sink node is in the first position of Eulerian cycle
    if tup.1 == cycle[0] {
        return cycle[1..].to_vec();
    }
    let mut path = cycle[idx + 1..cycle.len() - 1].to_vec();
    let rest = cycle[..idx + 1].to_vec();
    path.extend(rest);

    return path;
}
