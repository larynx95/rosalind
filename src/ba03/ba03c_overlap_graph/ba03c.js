/*
Rosalind: BA3C
Construct the Overlap Graph of a Collection of k-mers

In this chapter, we use the terms prefix and suffix
to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.

Given an arbitrary collection of k-mers Patterns,
we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern'
by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern').
The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

Overlap Graph Problem
Construct the overlap graph of a collection of k-mers.

Given:
A collection Patterns of k-mers.

Return:
The overlap graph Overlap(Patterns), in the form of an adjacency list.

Sample Dataset
ATGCG
GCATG
CATGC
AGGCA
GGCAT

Sample Output (1:1)
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * HERE: Construct OverapGraph (BA3C)
      -> NEXT: Reconstruct a String from its k-mer Composition (BA3H)

Plan 1.
  - deep copy patterns, create two list
  - use nested loop to make adjacent list

═════════════════════════════════════════════════

References:
- two forms of Graph: (a) adjacent list, (b) adjacent matrix
*/

// BA03C: overlap graph
function overlap_graph(texts) {
    /**
     * [str] -> {str:[str]}
     * returns a result of OverlapGraph(Patterns)
     */
    let dic = {};
    for (const source of texts) {
        for (const target of texts) {
            if (source.slice(1) === target.slice(0, target.length - 1)) {
                if (!(source in dic)) dic[source] = [target];
                else dic[source].push(target);
            }
        }
    }
    return dic;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03c.txt').toString().split("\n");

    const startTime = performance.now()
    const graph = overlap_graph(lines);
    for (const key in graph) {
        let line = "";
        line += key + " -> ";
        for (const val of graph[key]) {
            line += val + " ";
        }
        console.log(line);
    }
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()