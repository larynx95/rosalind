/*
Rosalind: BA3K
Generate "Contigs" from a Collection of Reads

Even after read breaking, most assemblies still have gaps in k-mer coverage,
causing the de Bruijn graph to have missing edges, and so the search for an Eulerian path fails.
In this case, biologists often settle on assembling "contigs" (long, contiguous segments of the genome)
rather than entire chromosomes.

For example, a typical bacterial sequencing project may result in about a hundred contigs,
ranging in length from a few thousand to a few hundred thousand nucleotides.
For most genomes, the order of these contigs along the genome remains unknown.
Needless to say, biologists would prefer to have the entire genomic sequence,
but the cost of ordering the contigs into a final assembly
and closing the gaps using more expensive experimental methods is often prohibitive.

Fortunately, we can derive contigs from the "de Bruijn" graph.
A path in a graph is called "non-branching" if "in(v) = out(v) = 1" for each intermediate node v of this path,
i.e., for each node except possibly the starting and ending node of a path.
A "maximal non-branching path" is a non-branching path that cannot be extended into a longer non-branching path.
We are interested in these paths because the strings of nucleotides
that they spell out must be present in any assembly with a given k-mer composition.
For this reason,
contigs correspond to strings spelled by "maximal non-branching paths" in the "de Bruijn graph".

Contig Generation Problem
Generate the contigs from a collection of reads (with imperfect coverage).

Given: A collection of k-mers Patterns.

Return: All contigs in DeBruijn(Patterns). (You may return the strings in any order.)

Sample Dataset
ATG
ATG
TGT
TGG
CAT
GGA
GAT
AGA

Sample Output
AGA ATG ATG CAT GAT TGGA TGT

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
    * PREV: Read-Pair (BA3J)
                                         "
              not every Eulerian path in the paired de Bruijn graph
         constructed from a (k, d)-mer composition spells out a solution of
                 the String Reconstruction from Read-Pairs Problem.
                                         "
      -> Construct a String Spelled by a Gapped Genome Path (BA3L)
      -> Generate All Maximal Non-Branching Paths in a Graph (BA3M)
      ↓
    * HERE: Generate "Contigs" from a Collection of Reads (BA3K)

Info.
  ╔═════════════════════════════════════════════════════════════════════════════╗
  ║ MAXIMALNONBRANCHINGPATHS(Graph)                                             ║
  ║   Paths <- empty list                                                       ║
  ║   for each node v in Graph                                                  ║
  ║     if v is not a 1-in-1-out node                                           ║
  ║       if OUT(v) > 0                                                         ║
  ║         for each outgoing edge (v, w) from v                                ║
  ║           NonBranchingPath <- the path consisting of the single edge (v, w) ║
  ║           while w is a 1-in-1-out node                                      ║
  ║             extend NonBranchingPath by the outgoing edge (w, u) from w      ║
  ║             w <- u                                                          ║
  ║           add NonBranchingPath to the set Paths                             ║
  ║   for each isolated cycle Cycle in Graph                                    ║
  ║     add Cycle to Paths                                                      ║
  ║   return Paths                                                              ║
  ╚═════════════════════════════════════════════════════════════════════════════╝

  * sample dataset
    ATG  AT->TG    AT:[TG,TG]    CA->AT=>TG->GT
    ATG  AT->TG                      ↑   ↓
    TGT  TG->GT    TG:[GT,GG]    AG->GA<-GG      AT->TG
    TGG  TG->GG                                  CA->AT  AT->TG  TG->GT
    CAT  CA->AT    CA:[AT]                       AG->GA  GA<-AT  TG->GG->GA
    GGA  GG->GA    GG:[GA]
    GAT  GA->AT    GA:[AT]
    AGA  AG->GA    AG:[GA]
    contigs: CAT, AGA, ATG, ATG, GAT, TGT, TGGA

  * another example
    DeBruijn_3(TAATGCCATGGGATGTT)                                              CC
    3-mers:                                                                   ↙ ↖
    TAA TA->AA    TA: [AA]                     CC                            CA   GC
    AAT AA->AT    AA: [AT]                   ↙  ↖                          ↓    ↑
    ATG AT->TG    AT: [TG,TG,TG]            CA    GC                         AT   TG
    TGC TG->GC    TG: [GC,GG,GT]            ↓     ↑                          AT->TG
    GCC GC->CC    GC: [CC]          TA->AA->AT ⇛ TG->GT->TT     TA->AA->AT  AT->TG   TG->GT->TT
    CCA CC->CA    CC: [CA]                  ↑     ↓                          AT->TG
    CAT CA->AT    CA: [AT]                  GT <- GG↺
    ATG AT->TG    GG: [GG,GA]                                   AT<-GA<-GG  TG->GG   GG↺
    TGG TG->GG    GA: [AT]
    GGG GG->GG    GT: [TT]
    GGA GG->GA
    GAT GA->AT
    ATG AT->TG
    TGT TG->GT
    GTT GT->TT
    contigs: TAAT, TGTT, TGCCAT, ATG, ATG, ATG, TGG, GGG, GGAT

═════════════════════════════════════════════════

╔═══════════════════════════════════════════════════════════╗
║ def allpaths(graph, start, sink, path=[]):                ║
║     """                                                   ║
║     ({a:[a]},a,a,[a]) -> [[a]]                            ║
║     returns all path from 'start' node to 'sink' node     ║
║     """                                                   ║
║     path = path + [start]                                 ║
║     if start not in graph:                                ║
║         return [path]                                     ║
║     paths = []                                            ║
║     for node in graph[start]:                             ║
║         if node not in path:                              ║
║             newpaths = allpaths(graph, node, sink, path)  ║
║             for newpath in newpaths:                      ║
║                 if newpath[-1] == sink:                   ║
║                     paths.append(newpath)                 ║
║     return paths                                          ║
╚═══════════════════════════════════════════════════════════╝
*/
// #!/usr/bin/env javascript

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

// BA03E: De Bruijn graph, from k-mers
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

// BA03M: maximum non-branching path
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
        const nodes = Object.keys(graph).filter(e => graph[e].length !== 0);  // some key has [] as its value
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

// main
function main() {
    const fs = require('fs');
    const kmers = fs.readFileSync('/home/wsl/rosalind/data/ba03k.txt').toString().split("\n");

    const startTime = performance.now()
    const graph = debruijn_from_patterns(kmers);
    const paths = max_nonbrancing_path(graph);
    console.log(paths.join(' '));
    console.log(`${performance.now() - startTime} milliseconds`)

    /*
    var stream = fs.createWriteStream("../data/ba3k_output.txt");
    stream.once('open', function (fd) {
        for (const contig of paths) {
            stream.write(contig);
            stream.write('\n');
        }
        stream.end();
    });
    */
}

// execute main function
main()