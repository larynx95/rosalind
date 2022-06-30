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
    CTTA -> CTT,TTA
    ACCA -> ACC,CCA
    TACC -> TAC,ACC
    GGCT -> GGC,GCT
    GCTT -> GCT,CTT
    TTAC -> TTA,TAC
    prefixes: [CCT,ACC,TAC,GGC,GCT,TTA]
    suffixes: [TTA,CCA,ACC,GCT,CTT,TAC]

Plan 2.
  * using Eulerian path (BA3G) algorithm
  * multiple strings -> De Bruijn graph -> Eulerian path

═════════════════════════════════════════════════

References:
-
*/

// BA03E: De Bruijn graph, from k-mers
function debruijn_from_kmers(kmers) {
    /**
     * [str] -> {str:[str]}
     * returns a De Bruijn graph
     */
    let graph = {};
    for (const kmer of kmers) {
        const prefix = kmer.slice(0, kmer.length - 1);
        const suffix = kmer.slice(1);
        if (!(prefix in graph)) graph[prefix] = [suffix];
        else graph[prefix].push(suffix);
    }
    return graph;
}

// simple walk
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

// BA03F: Eulerian cycle
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

// BA01G: Eulerian path
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

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03h.txt').toString().split("\n");
    const k = parseInt(lines[0]);
    const kmers = lines.slice(1);

    const startTime = performance.now()
    const graph = debruijn_from_kmers(kmers);
    const path = eulerian_path(graph);
    let answer = path[0];
    for (const elem of path.slice(1)) {
        answer += elem[elem.length - 1];
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()