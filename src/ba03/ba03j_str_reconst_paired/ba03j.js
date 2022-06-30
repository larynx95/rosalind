/*
Rosalind: BA3J
Reconstruct a String from its Paired Composition

Since increasing read length presents a difficult experimental problem,
biologists have suggested an indirect way of increasing read length by generating read-pairs,
which are pairs of reads separated by a fixed distance d in the genome.

You can think about a read-pair as a long "gapped" read of length k + d + k
whose first and last k-mers are known but whose middle segment of length d is unknown.
Nevertheless, read-pairs contain more information than k-mers alone,
and so we should be able to use them to improve our assemblies.
If only you could infer the nucleotides in the middle segment of a read-pair,
you would immediately increase the read length from k to 2 · k + d.

Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d.
We use the notation (Pattern1|Pattern2) to refer to a a (k,d)-mer whose k-mers are Pattern1 and Pattern2.
The (k,d)-mer composition of Text, denoted PairedCompositionk,d(Text),
is the collection of all (k,d)- mers in Text (including repeated (k,d)-mers).

String Reconstruction from Read-Pairs Problem
Reconstruct a string from its paired composition.

Given: Integers k and d followed by a collection of paired k-mers PairedReads.

Return: A string Text with (k, d)-mer composition equal to PairedReads.
        (If multiple answers exist, you may return any one.)

Sample Dataset
4 2
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT

Sample Output
GTGGTCGTGAGATGTTGA

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
    * HERE: Read-Pair (BA3J)
                                         "
              not every Eulerian path in the paired de Bruijn graph
         constructed from a (k, d)-mer composition spells out a solution of
                 the String Reconstruction from Read-Pairs Problem.
                                         "
      -> Construct a String Spelled by a Gapped Genome Path (BA3L)
      -> PREV: Generate All Maximal Non-Branching Paths in a Graph (BA3M)
      ↓
    * NEXT: Generate "Contigs" from a Collection of Reads (BA3K)

Info.
  * Terminology
    - Prefix, Suffix
      ; Suffix((TAA|GCC)) = Prefix((AAT|CCA)) = (AA|CC)
      ; Prefix((ATG|CAT)) = Suffix((AAT|CCA)) = (AT|CA)

    - (k,d)-mer
      ; 'k-mer' is a pattern with length k
      ; 'd' is distance between Pattern1 and Pattern2 in (Pattern1 | Pattern2)
      ; e.g. TAA-GCC in TAATGCCATGGGATGTT

    - (Pattern1 | Patern2)
      ; (k,d)-mer whose k-mers are pattern1 and Pattern2

    - PairedComposition_(k,d)(Text)
      ; a collection PairedReads - all (k,d)-mers in Text
      ; e.g. PairedComposition_(3,1)(TAATGCCATGGGATGTT)
        (AAT|CCA) (ATG|CAT) (ATG|GAT) (CAT|GGA) (CCA|GGG) ...

    - CompositionGraph_(k,d)(Text)
      ; graph consisting of |Text| - (2k+d) +1 isolated edges
        that are labeled by the (k,d)-mers in Text,
        and whose nodes are labeled by the prefixes and suffixes of these labels
      ; (Important) There may be overlapping nodes.

    - PathGraph_(k,d)(Text)
      ; path formed by |Text| - (2k+d) + 1 edges
      ; e.g. PathGraph_(3,1)(TAATGCCATGGGATGTT)
           TAA    AAT    ATG
           GCC    CCA    CAT
        TA ->  AA ->  AT ->  ...
        GC     CC     CA

  * example - 'PairedComposition_(3,1)(Text)'
    index: 01234567890123456       length of a pair: k + d + k = 2k + d = 2*3+1 = 7
           TAATGCCATGGGATGTT       number of pairs : |Text| - (2k+d) + 1

  * example
    TAA GCC
     AAT CCA            (k,d)-mers
      ATG CAT           Suffix((TAA|GCC)) = Prefix((AAT|CCA)) = (AA|CC)
       TGC ATG                   ^^  ^^             ^^  ^^
        GCC TGG         Prefix((ATG|CAT)) = Suffix((AAT|CCA)) = (AT|CA)
         CCA GGG                ^^  ^^               ^^  ^^
          CAT GGA
           ATG GAT
            TGG ATG
             GGG TGT
              GGA GTT
    TAATGCCATGGGATGTT

  * split example
    TAA GCC ┌ Prefix: (TA, GC)
            └ Suffix: (AA, CC)
    AAT CCA ┌ Prefix: (AA, CC)
            └ Suffix: (AT, CA)
    ATG CAT ┌ Prefix: (AT, CA)
            └ Suffix: (TG, AT)
    TGC ATG ┌ Prefix: (TG, AT)
            └ Suffix: (GC, TG)
    GCC TGG ┌ Prefix: (GC, TG)
            └ Suffix: (CC, GG)
    CCA GGG ┌ Prefix: (CC, GG)
            └ Suffix: (CA, GG)
    CAT GGA ┌ Prefix: (CA, GG)
            └ Suffix: (AT, GA)
    ATG GAT ┌ Prefix: (AT, GA)
            └ Suffix: (TG, AT)
    TGG ATG ┌ Prefix: (TG, AT)
            └ Suffix: (GG, TG)
    GGG TGT ┌ Prefix: (GG, TG)
            └ Suffix: (GG, GT)
    GGA GTT ┌ Prefix: (GG, GT)
            └ Suffix: (GA, TT)
    TAATGCCATGGGATGTT

  * How to gluing PathGraph_(k,d)(Text)
    4 2          prefix     suffix        Eulerian path
    --------------------------------------------------------
    GTGG|GTGA -> (GTG|GTG)  (TGG|TGA)     GTG--GTG
    TGGT|TGAG -> (TGG|TGA)  (GGT|GAG)      TGG--TGA
    GGTC|GAGA -> (GGT|GAG)  (GTC|AGA)       GGT--GAG
    GTCG|AGAT -> (GTC|AGA)  (TCG|GAT)        GTC--AGA
    TCGT|GATG -> (TCG|GAT)  (CGT|ATG)         TCG--GAT
    CGTG|ATGT -> (CGT|ATG)  (GTG|TGT)          CGT--ATG
    GTGA|TGTT -> (GTG|TGT)  (TGA|GTT)           GTG--TGT
    TGAG|GTTG -> (TGA|GTT)  (GAG|TTG)            TGA--GTT
    GAGA|TTGA -> (GAG|TTG)  (AGA|TGA)             GAG--TTG
                                                   AGA--TGA
                                          GTGGTCGTGAGA
                                          4+2   GTGAGATGTTGA
                                                012345678901234567
  * TODO:
    "... not every Eulerian path in the paired de Bruijn graph
     constructed from a (k,d)-mer composition
     spells out a soultion of the String Reconstruction from Read-Pairs problem..."

Plan 1.
  * steps:
    - get 'PairedDeBruijnGraph_(k,d)(Text)'
    - find end node and start node, connect two node to make Eulerian path to Eulerian Cycle
    - get 'EulerianCycle_(k,d)(Text)'
    - disconnect end and start nodes, get 'EulerianPath_(k,d)(Text)' --> 'PathGraph_(k,d)(Text)'
    - glue 'PathGraph_(k,d)(Text)' to get string
      string1 + string2[k+d:]
      => Failed! TODO: Why?

  * (CAUTION) This approach is not always correct!
    "every Eulerian path in the de Bruijn graph constructed from a k-mer composition
     spells out a solution of the String Reconstruction Problem.
     But is this the case for the paired de Bruijn graph?"

    "not every Eulerian path in the paired de Bruijn graph constructed
     from a (k, d)-mer composition spells out a solution of
     the String Reconstruction from Read-Pairs Problem."

  * figure 3.36

            AG    A    CA
       AG   TG ↙ T ↖ CT   CT
       AG   ↙        ↖    CA
     A -> G  ⇛  GC ⇛    C  ->  T
     A    G      GC       C      A
            ↖         ↙
          TG  ↖  T  ↙  CT
          TG      T      CT

  (1) Eulerian Path 1
    (AG|AG) -> (GC|GC) -> (CT|CT) -> (TG|TG) -> (GC|GC) -> (CA|CT) -> (AG|TG) -> (GC|GC) -> (CT|CA)
    (AG|AG)
     (GC|GC)
      (CT|CT)
       (TG|TG)
        (GC|GC)
         (CA|CT)
          (AG|TG)
           (GC|GC)
            (CT|CA)
     AGCTGCAGCT
           AGCTGCTGCA    <-- The overlapping part depends on the path.

  (2) Eulerian Path 2
    (AG|AG) -> (GC|GC) -> (CA|CT) -> (AG|TG) -> (GC|GC) -> (CT|CT) -> (TG|TG) -> (GC|GC) -> (CT|CA)
    (AG|AG)
     (GC|GC)
      (CA|CT)
       (AG|TG)
        (GC|GC)
         (CT|CT)
          (TG|TG)
           (GC|GC)
            (CT|CA)
     AGCAGCTGCT
        AGCTGCTGCA     <-- The overlapping part depends on the path.

Plan 2.
  * find the index where string overapping starts.
    then string = string1[:index] + string2
  * finding the index or indices
    - built-in 'find()'
    - using BA1D: Find All Occurrences of a Pattern in a String

═════════════════════════════════════════════════

References:
- add vs update in set operations in python
  https://stackoverflow.com/questions/28845284/add-vs-update-in-set-operations-in-python
- How to find all occurrences of a substring?
  https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
  "Python is "battery-included language"
*/

// construct De Bruijn graph, from string lines
function debruijn_paired(lines) {
    /**
     * [str] -> {str:[str]}
     * returns a paired De Bruijn graph
     * In JavaScript, object' key is only string type
     *                non-primitive type key is not unique in Map/Set  --> This is annoying.
     */
    let graph = {};  // Object's property is string type, not else. Annoying. I miss dictionary of Python.
    for (const line of lines) {
        if (line == "") break;  // if the last line is empty space
        const [fst, snd] = line.split('|');
        const prefix = fst.slice(0, fst.length - 1) + "|" + snd.slice(0, snd.length - 1);
        const suffix = fst.slice(1) + "|" + snd.slice(1);
        if (!(prefix in graph)) graph[prefix] = [suffix];
        else graph[prefix].push(suffix);
    }
    return graph;
}

// simple walk, construct graph
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
     * returns a Eluerian cycle, destructive (BA3F)
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
     * returns In/Out-degrees (BA3F)
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

// BA03G: Eulerain path
function eulerian_path(graph) {
    /**
     * {str:[str]} -> [str]
     * returns a Eulerian path (BA3G)
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

// BA03J: string reconstruction from paired compositions
function str_reconst_paired_composition(graph, k, d) {
    /**
     * [str] -> str
     * returns a reconstructed string from a list of strings
     */
    const path = eulerian_path(graph);
    let fst = path[0].slice(0, k - 1);
    let snd = path[0].slice(path[0].length - (k - 1));
    for (const pair of path.slice(1)) {
        fst += pair[k - 2];
        snd += pair[pair.length - 1];
    }
    return fst + snd.slice(k + d);
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03j.txt').toString().split("\n");
    const [k, d] = lines[0].split(/\s/).map(Number);
    const graph = debruijn_paired(lines.slice(1));

    const startTime = performance.now();
    console.log(str_reconst_paired_composition(graph, k, d));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
