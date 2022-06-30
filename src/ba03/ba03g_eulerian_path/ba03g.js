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

Plan 2.
  * find start, sink nodes
  * walk from strat node to sink node

═════════════════════════════════════════════════

References:
-
*/

// simple walk, constructing graph
function walk(graph, node) {
    /**
     * ({str:[str]},str) -> [str]
     * returns a path or cycle by simple random walk
     */
    let cycle = [node];
    while (graph[node].length != 0) {
        const values = graph[node];
        const target = values[Math.floor(Math.random() * values.length)];
        cycle.push(target);
        graph[node] = values.filter(e => e != target);
        node = target;
    }
    return cycle;
}

// BA03F: Hierholzer's algorithm
function eulerian_hierholzer(graph) {
    /**
     * {a:[a]} -> [a]
     * returns a Eluerian cycle, destructive
     * >>> var graph = {'0':['3'],'1':['0'],'2':['1','6'],'3':['2'],'4':['2'],'5':['4'],'6':['5','8'],'7':['9'],'8':['7'],'9':['6']};
     * >>> eulerian_hierholzer(graph)
     *     ['6','5','4','2','1','0','3','2','6','8','7','9','6']
     */
    // first, get a cycle by random walk
    let nodes = Object.keys(graph);
    let source = nodes[Math.floor(Math.random() * nodes.length)];  // initial source node, randomly selected
    let cycle = walk(graph, source);
    // if there's remaining unvisited nodes, get another cycle, merge it to the previous cycle, again and again
    while (Object.keys(graph).map(e => graph[e].length === 0).every(Boolean) != true) {
        source = cycle.filter(e => graph[e].length != 0)[0];
        const new_cycle = walk(graph, source);
        const idx = cycle.indexOf(source);
        cycle = [...cycle.slice(0, idx + 1), ...new_cycle.slice(1), ...cycle.slice(idx + 1)];
    }
    return cycle;
}

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

// BA03G: Eulerian path, by disconnecting Eulerian cycle
function eulerian_path(graph) {
    /**
     * {str:[str]} -> [str]
     * returns a Eulerian path
     */
    // get In/Out-degree dictionary (node: [In, Out])
    let iodic = {};
    for (const key in graph) {
        // count In-degrees
        for (const val of graph[key]) {
            if (!(val in iodic)) iodic[val] = [1, 0];
            else iodic[val][0]++;
        }
        // count Out-degrees
        if (!(key in iodic)) iodic[key] = [0, graph[key].length];
        else iodic[key][1] += graph[key].length;
    }
    // find start, sink nodes
    let start, sink;
    for (const key in iodic) {
        if (iodic[key][0] - iodic[key][1] < 0) start = key;
        else if (iodic[key][0] - iodic[key][1] > 0) sink = key;
    }
    // connect sink to start: create a new edge between start and sink nodes
    if (!(sink in graph)) graph[sink] = [start];
    else graph[sink].push(start);
    // create a Eulerian cycle
    let cycle = eulerian_hierholzer(graph);
    cycle = cycle.slice(1);
    // rotate, disconnect
    while (cycle[cycle.length - 1] !== sink) {
        cycle.push(cycle[0]);
        cycle = cycle.slice(1);
    }
    return cycle;
}

// helper fx
function read_graph(lines) {
    /**
     * [str] -> {int:[int]}  <-- This is wrong in Javascript! Object key must be a string type.
     * [str] -> {str:[str]}  <-- This is correct. Very strange.
     * returns a Eulerian graph ("node: [In-degrees, Out-degrees]"")
     */
    let graph = {};
    for (const line of lines) {
        const [source, val_str] = line.split(' -> ');
        let vals = [];
        for (const val of val_str.split('\,')) {
            vals.push(val);
        }
        graph[source] = vals;
    }
    return graph;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03g.txt').toString().split("\n");
    let graph = read_graph(lines);

    const startTime = performance.now()
    const path = eulerian_path(graph);
    console.log(path.join('->'));
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()