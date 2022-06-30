/*
Rosalind: BA3F (difficulty: 4/5)
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

╔══════════════════════════════════════════════════════════════════════╗
║ EULERIANCYCLE(Graph)                                                 ║
║     form a cycle Cycle by randomly walking in Graph                  ║
║         (don't visit the same edge twice!)                           ║
║     while there are unexplored edges in Graph                        ║
║         select a node newStart in Cycle with still unexplored edges  ║
║         form Cycle' by traversing Cycle (starting at newStart)       ║
║           and then randomly walking                                  ║
║         Cycle <- Cycle'                                              ║
║     return Cycle                                                     ║
╚══════════════════════════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║
║         v <- a node with indegree larger than 1 in G        ║
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║
║                 if NewGraph is connected                    ║
║                     add NewGraph to ALlGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝

Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An "Eulerian directed graph", in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3
1 -> 0
2 -> 1,6
3 -> 2
4 -> 2
5 -> 4
6 -> 5,8
7 -> 9
8 -> 7
9 -> 6

Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

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
      -> PREV: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> HERE: Eulerian Cycle (BA3F)
      -> NEXT: Eulerian Path (BA3G)

Info.
1. Hierholzer's algorithm
  * steps:
    ; just walk through the graph, and get a circle1
    ; if the circle1 is complete whole large cycle without remaining nodes in graph, that's Eulerian cycle
    ; if there're remaining unvisited node in graph yet, just do another walking, get another cycle2
      but in this time, start node should be a node in previous cycle, and it has one or more out-degrees
    ; patch cycle2 to cycle1, repeat previous steps
              4 < 5             4 < 5             4 * 5             4 * 5
              v   ^             v   ^             *   *             *   *
          1 < 2 > 6 > 8     1 *[2]> 6 > 8     1 * 2 * 6 > 8     1 * 2 *[6]* 8
          v   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
          0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
                            10321             10321             103265421
                                                 26542              68796
                                              103265421         1032687965421

2. Rust specific things
  * function arguments:
    - reference vs. value
    - &String vs. &str
  * life time
    - variable declared in local scope: type, life time annotation?
  * modifying value in hashMap<String, Vec<String>>

═════════════════════════════════════════════════

References:
- How can I update a value in a mutable HashMap?
  https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap
- if argument is '&mut HashMap', what is the type of 'map.get(&key).unwrap()'?
*/

use rand::Rng;
use std::collections::HashMap;
use std::time::Instant;
use std::{
    fs::File,
    io::{prelude::*, BufReader},
};

// walk
// To be considered ...
// 1. type of arguments
//   - Function 'walk' changes it's arguments repeatedly, so mutable references were used here.
//   - Original graph should be modified, so mutable reference '&mut' is used.
// 2. considering life-time
//   - 'source', 'target' is local variable in while loop. So I used '&String' not '&str'.
// 3. operations
//   - HashMap iteration
//   - get value by key mutably
//   - 'clone()' to use moved variable
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
        println!("{:?}", graph);
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

// main
fn main() {
    let lines = lines_from_file("/home/wsl/rosalind/data/ba03f.txt");

    // check time
    let start = Instant::now();
    let mut graph: HashMap<String, Vec<String>> = read_graph(&lines);
    let cycle = eulerian_cycle(&mut graph);
    println!("{}", cycle.join("->"));
    println!("Successful! {:?}", start.elapsed());
}

// helper fx: read a text file, create a Vector of strings
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
