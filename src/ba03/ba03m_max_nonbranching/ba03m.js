/*
Rosalind: BA3M
Generate All Maximal Non-Branching Paths in a Graph

A node v in a directed graph Graph is called a 1-in-1-out node
if its indegree and outdegree are both equal to 1, i.e., in(v) = out(v) = 1.
We can rephrase the definition of a "maximal non-branching path"
from the main text as a path whose internal nodes are 1-in-1-out nodes
and whose initial and final nodes are not 1-in-1-out nodes.
Also, note that the definition from the main text does not handle the special case
when Graph has a connected component that is an isolated cycle,
in which all nodes are 1-in-1-out nodes.

The MaximalNonBranchingPaths pseudocode below generates all non-branching paths in a graph.
It iterates through all nodes of the graph that are not 1-in-1-out nodes
and generates all non-branching paths starting at each such node.
In a final step, MaximalNonBranchingPaths finds all isolated cycles in the graph.

╔═══════════════════════════════════════════════════════════════════════════════╗
║  MaximalNonBranchingPaths(Graph)                                              ║
║    Paths <- empty list                                                        ║
║    for each node v in Graph                                                   ║
║      if v is not a 1-in-1-out node                                            ║
║        if out(v) > 0                                                          ║
║          for each outgoing edge (v, w) from v                                 ║
║            NonBranchingPath <- the path consisting of the single edge (v, w)  ║
║              while w is a 1-in-1-out node                                     ║
║                extend NonBranchingPath by the outgoing edge (w, u) from w     ║
║                  w <- u                                                       ║
║              add NonBranchingPath to the set Paths                            ║
║    for each isolated cycle Cycle in Graph                                     ║
║      add Cycle to Paths                                                       ║
║    return Paths                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝

Maximal Non-Branching Path Problem
Find all maximal non-branching paths in a graph.

Given: The adjacency list of a graph whose nodes are integers.

Return: The collection of all maximal non-branching paths in the graph.

Sample Dataset
1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6

Sample Output
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7

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
      -> Eulerian Cycle (BA3F)
        -> k-universal string problem (BA3I)
      -> Eulerian Path (BA3G)
      ↓
    * NEXT: Read-Pair (BA3J)
                                         "
              not every Eulerian path in the paired de Bruijn graph
         constructed from a (k, d)-mer composition spells out a solution of
                 the String Reconstruction from Read-Pairs Problem.
                                         "
      -> PREV: Construct a String Spelled by a Gapped Genome Path (BA3L)
      -> HERE: Generate All Maximal Non-Branching Paths in a Graph (BA3M)

Plan 1.
  * sample dataset
    [1] -> 2 -> [3] -> [4]      6 -> 7 -> 6
                 ↓
                [5]
    v           v       v    <-- not 1-in-1-out nodes: 1,3,4,5
    1:[2] 2:[3] 3:[4,5] 6:[7] 7:[6]

  * non-branching paths
    1 -> 2 -> 3 ┌ 3 -> 4
                └ 3 -> 5
    6 -> 7 -> 6

  * degrees
    {1:[2], 2:[3], 3:[4,5], 6:[7], 7:[6]}            <-- graph
    {1:[0,1], 2:[0,1], 3:[0,2], 6:[0,1], 7:[0,1]}    <-- out-degrees: {source:[-, len(graph[source])]}
    {1:[0,1], 2:[0,1], 3:[1,2], 4:[1,0], 5:[1,0], 6:[1,1], 7:[1,1]}  <-- in-degrees
*/

// get degrees
function get_degrees(graph) {
    /**
     * {str:[a]} -> {str:[int,int]}       // no tuple type in Javascript, object property is always string
     * returns In/Out-degrees
     * >>> get_degrees({ 0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6] })
     *     {'0':[1,1],'1':[1,1],'2':[2,2],'3':[1,1],'4':[1,1],'5':[1,1],'6':[2,2],'7':[1,1],'8':[1,1],'9':[1,1]}
     */
    let dic = {};
    for (const key in graph) {
        // count In-degrees
        for (const val of graph[key]) {
            if (!(val in dic)) dic[val] = [1, 0];
            else dic[val][0]++;
        }
        // count Out-degrees
        if (!(key in dic)) dic[key] = [0, graph[key].length];
        else dic[key][1] += graph[key].length;
    }
    return dic;
}

// BA03M: maximum non-branching path (contig)
function max_nonbrancing_path(graph) {
    /**
     * {str:[str]} -> [str]
     * returns a list of contigs (BA3K)
     */
    // I. get a dictionary of degrees
    const degrees = get_degrees(graph);
    // II. get non-breaking paths (contigs)
    let paths = [];
    // A. linear paths
    for (const key in graph) {
        if (!(degrees[key][0] == 1 && degrees[key][1] == 1) && degrees[key][1] > 0) {
            for (const val of graph[key]) {
                let nbpath = key + val[val.length - 1];  // nbpath = [key, val];
                let cur = val;
                graph[key] = graph[key].filter(e => e !== val);
                while (degrees[cur][0] == 1 && degrees[cur][1] == 1) {
                    const target = graph[cur].pop();
                    nbpath += target[target.length - 1];  // nbpath.push(target);
                    cur = target;
                }
                paths.push(nbpath);
            }
        }
    }
    // B. cycles
    while (Object.keys(graph).map(e => graph[e].length === 0).every(Boolean) != true) {
        const nodes = Object.keys(graph).filter(e => graph[e].length !== 0);
        const start = nodes[Math.floor(Math.random() * nodes.length)];
        let target = graph[start].pop();
        let cycle = start + target[target.length - 1];
        let cur = target;
        while (cur !== start) {
            target = graph[cur].pop();
            cycle += target[target.length - 1];
            cur = target;
        }
        paths.push(cycle);
    }
    return paths;
}

// helper fx
function read_graph(lines) {
    /**
     * [str] -> {int:[int]}  <-- This is wrong in Javascript! Object key must be a string type.
     * [str] -> {str:[str]}  <-- This is correct. Very strange.
     * returns a Eulerian graph
     */
    let graph = {};
    for (const line of lines) {
        const [source, val_str] = line.split(' -> ');
        let vals = [];
        for (const val of val_str.split(',')) {
            vals.push(val);
        }
        graph[parseInt(source)] = vals;
    }
    return graph;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03m.txt').toString().split("\n");

    const startTime = performance.now()
    const graph = read_graph(lines);
    const paths = max_nonbrancing_path(graph);
    for (const path of paths) {
        console.log(Array.from(path).join(' -> '));
    }
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()