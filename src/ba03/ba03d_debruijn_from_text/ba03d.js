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
TCT -> CTA,CTC
TTC -> TCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> PREV: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> HERE: Construct De Bruijn Graph (BA3D)
      -> NEXT: Construct De Bruijn Graph with k-mers (BA3E)

Info.
  * Example (Overlap graph):
    TAATGCCATGGGATGTT
    TAA AAT ATG TGC GCC CCA CAT ATG TGG GGG GGA GAT ATG TGT GTT

    Example (De Bruijn graph):
    edges:  TAA AAT ATG TGC GCC CCA CAT ATG TGG GGG GGA GAT ATG TGT GTT
    nodes: TA->AA->AT->TG->GC->CC->CA->AT->TG->GG->GG->GA->AT->TG->GT->TT
    repeat:         *   ^               *   ^   v   v       *   ^

Plan 1.
- AAGATTCTCTAC
  => k-mers    :  AAGA    AGAT    GATT    ATTC    TTCT    TCTC    CTCT    TCTA    CTAC
                  /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \
  => (k-1)-mers: AAG AGA AGA GAT GAT ATT ATT TTC TTC TCT TCT CTC CTC TCT TCT CTA CTA TAC
                 1   2   2   2   2   2   2   2   2   4   4   2   2   4   4   2   2   1
  => (AAG,AGA) (AGA,GAT) (GAT,ATT) (ATT,TTC) (TTC,TCT) (TCT,CTC) (CTC,TCT) (TCT,CTA) (CTA,TAC)
  (AAG,AGA)
      (AGA,GAT)
          (GAT,ATT)
              (ATT,TTC)
                  (TTC,TCT)
                      (TCT,CTC)
                          (CTC,TCT)
                      (TCT,CTA)       <-- TCT again!
                          (CTA,TAC)

═════════════════════════════════════════════════

References:
- How To Pronounce "De Bruijn"?  [debroin]
  https://www.biostars.org/p/7186/
*/

// BA03A: all k-mer
function* all_kmers_gen(text, k) {
    /**
     * (str,int) -> generator
     * returns a generator of all k-mers
     */
    for (let i = 0; i < text.length; i++) yield text.slice(i, i + k);
}

// BA03A: all k-mer
function all_kmers(text, k) {
    /**
     * (str,int) -> [str]
     * returns a list of all k-mers
     */
    return [...Array(text.length - k + 1)].map((_, i) => e = text.slice(i, i + k));
}

// BA03D: De Bruijn graph, from a text
function debruijn_from_text(text, k) {
    /**
     * (int,str) -> {str:[str]}
     * returns a De Bruijn graph from a text (BA3D)
     */
    let graph = {};
    for (let i = 0; i < text.length - k + 1; i++) {
        const kmer = text.slice(i, i + k);
        const prefix = kmer.slice(0, kmer.length - 1);
        const suffix = kmer.slice(1);
        if (!(prefix in graph)) graph[prefix] = [suffix];
        else graph[prefix].push(suffix);
    }
    return graph;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03d.txt').toString().split("\n");
    const k = parseInt(lines[0]);
    const text = lines[1];

    // display on screen
    const startTime = performance.now()
    const debruijn_graph = debruijn_from_text(text, k);
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