/*
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node is called an Eulerian path.

Eulerian Path Problem
Find an Eulerian path in a graph.

Given: A directed graph that contains an Eulerian path,
       where the graph is given in the form of an adjacency list.

Return: An Eulerian path in this graph.

Sample Dataset
0 -> 2
1 -> 3
2 -> 1
3 -> 0,4
6 -> 3,7
7 -> 8
8 -> 9
9 -> 6

Sample Output
6->7->8->9->6->3->0->2->1->3->4

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
      -> Construct De Bruijn Graph (BA3D)
      -> Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> PREV: Eulerian Cycle (BA3F)
        -> NEXT: k-universal string problem (BA3I)
      -> HERE: Eulerian Path (BA3G)

Info.
  * writing down graph to figure out what this exercise is:
    {0:[2], 1:[3], 2:[1], 3:[0, 4], 6:[3,7], 7:[8], 8:[9], 9:[6]}
                  4
                  ↑
          9 → 6 → 3 → 0       Not circular!!!
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
  * problems:
    - This graph is not circular graph.
      Can any node be a starting point (source) or an ending point (target)?
    - How can I find the end node?
      a. comparing keys and values in dictionary (De Bruijn graph)
         dict.keys:   0123  6789
         dict.values: 012334 789   --> 4 is the end!
      b. difference bw in-degree and out-degree
         In-degree must be larger than out-degree in the last node.
         0:0, 1:0, 2:0, 3:0, 4:1, 6:-1, 7:0, 8:0, 9:0  --> The node 4 is the last!

Plan 1.
  * re-using "EulerianCycle" function
    - Adding an extra edge from w (end node) to v (start node)
      transforms the Eulerian path into an Eulerian cycle.
    - convert Eulerian path to Eulerian cycle, then re-convert cycle to path

  * out-degree +1, in-degree -1
                  4
                  ↑
          9 → 6 → 3 → 0
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
    unique nodes: 0    1    2    3    4    6    7    8    9
    --------------------------------------------------------
      0:[2]       +1       -1
      1:[3]           +1        -1
      2:[1]           -1   +1
      3:[0,4]     -1            +2   -1
      6:[3,7]                   -1        +2   -1
      7:[8]                                    +1   -1
      8:[9]                                         +1   -1
      9:[6]                               -1             +1
    --------------------------------------------------------
                  0    0    0    0   -1   +1    0    0    0   --> start:6 end (sink):4
  * add edge (end, start) to graph, get Eulerian cycle: 4 -> 6
  * trim the last repetitive elem in cycle, then rotate cycle to get Eulerian path

═════════════════════════════════════════════════

References:
- Rust by Example: Tuples
  https://www.cs.brandeis.edu/~cs146a/rust/rustbyexample-02-21-2015/tuples.html
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
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03g.txt");

    // check time
    let start = Instant::now();
    let mut graph = read_graph(&lines);
    let path = eulerian_path(&mut graph);
    println!("{}", path.join("->"));
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
fn print_graph(graph: HashMap<String, Vec<String>>) {
    for (key, values) in graph {
        print!("{key} -> ");
        println!("{}", values.join(","));
    }
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
            if !graph.get(node).unwrap().is_empty() {
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
