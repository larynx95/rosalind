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
- TODO: Is De Bruijn Graph is the same as Eulerian Graph?
*/

// BA03E: De Bruijn graph, from patterns
function debruijn_from_patterns(patterns) {
    /**
     * [str] -> {str:[str]}
     * returns a De Bruijn graph from patterns (BA3E)
     */
    let graph = {};
    for (const pattern of patterns) {
        const prefix = pattern.slice(0, pattern.length - 1);
        const suffix = pattern.slice(1);
        if (!(prefix in graph)) graph[prefix] = [suffix];
        else graph[prefix].push(suffix);
    }
    return graph;
}

// main
function main() {
    const fs = require('fs');
    const patterns = fs.readFileSync('/home/wsl/rosalind/data/ba03e.txt').toString().split("\n");

    // display on screen
    const startTime = performance.now()
    const debruijn_graph = debruijn_from_patterns(patterns);
    for (const key in debruijn_graph) {
        let answer = key + ' -> ' + debruijn_graph[key].join(',');
        console.log(answer);
    }
    console.log(`${performance.now() - startTime} milliseconds`)

    /*
    // write to file
    var stream = fs.createWriteStream("../data/ba3d_output.txt");
    stream.once('open', function (fd) {
        for (const key in debruijn_graph) {
            let answer = key + ' -> ' + debruijn_graph[key].join(',') + '\n';
            stream.write(answer);
        }
        stream.end();
    });
    */
}

// execute main function
main()