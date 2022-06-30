/*
Rosalind: BA3I
Find a k-Universal Circular String

A k-universal circular string is a circular string
that contains every possible k-mer constructed over a given alphabet.

k-Universal Circular String Problem
Find a k-universal circular binary string.

Given: An integer k.

Return: A k-universal circular string.
(If multiple answers exist, you may return any one.)

Sample Dataset
4

Sample Output
0000110010111101

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
        -> HERE: k-universal string problem (BA3I)
      -> PREV: Eulerian Path (BA3G)
      ↓
    * Read-Pair (BA3J)
      -> NEXT: Construct a String Spelled by a Gapped Genome Path (BA3L)

Info.
  * sample dataset if k=3, k-universal string is "01000111"
    01           010    .
     10          .100   .
      00         . 000  .
       00        .  001 .
        01       .   011.
         11      .    111 <--
          11     .     11-0
           10    .      1-01
            01
    0100011101   01000111-01
    if linear    if cycle

Plan 1.
  * product
    k=1 [0,1] 2^1  [0,1]
    k=2 [0,1] 2^2  [00,01,10,11]
    k=3 [0,1] 2^3  [000,001,010,100,011,101,110,111]

  * steps:
    i.   create a list k-mers - cartesian product with "01"
         All (not some) of these binary k-mers should be used to reconstruct k-universal string.
    ii.  get De Bruijn graph from the list (exclude edges to node itself)
         split each k-mer into prefix and suffix (k-1)-mers
         then construct De Bruijn graph represented as dictionary
    iii. get Eulerian cycle from De Bruijn graph
         TODO: Is the graph always Eulerian cycle? Always balanced?
               Strangely, the graph is always Eulerian. Why?
         TODO: If the graph is not Eulerian, what shall I do?
    iv.  find out the way how to reconstruct string from Eulerian cycle
         BA3B string reconstruction function is not for circular. Be careful!
         - temp_string = cycle[0] + last characters from cycle[1:]
         - final_string = temp_string[:-(k-1)]

═════════════════════════════════════════════════

References:
- Generate permutations of JavaScript array
  https://stackoverflow.com/questions/37579994/generate-permutations-of-javascript-array
*/

// helper fx. : get all product
function product_rec(arr, k) {
    /**
     * ([str],int) -> [str]
     * returns a list of product from fiven list of characters
     * >>> product_rec(['0','1'], 2)
     *     ['00','01','10','11]
     */
    if (k == 1) return arr;
    let res = [];
    for (const elem of arr) {
        for (const s of product_rec(arr, k - 1)) {
            res.push(elem + s);
        }
    }
    return res;
}

// helper fx. : get all product, generator version
function* product_rec_gen(arr, k) {
    /**
     * ([str],int) -> generator
     * returns a list of product from fiven list of characters, generator version
     * I asked stackoverflow about how to implement this function.
     * https://stackoverflow.com/questions/72195200/javascript-getting-all-possible-k-length-strings-from-given-list-of-characters
     */
    for (let i = 0; i < arr.length; i++) {
        if (k === 1) {
            yield [arr[i]];
        } else {
            const rest = product_rec_gen(arr, k - 1);
            for (let el of rest) {
                yield arr[i] + el;
            }
        }
    }
}

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

// check if graph is Eulerian
function check_eulerian_cycle(graph) {
    /**
     * {str:[a]} -> bool
     * returns true if graph is Eulerian cycle, else false (all nodes with even degrees)
     * >>> check_eulerian_cycle({0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]})
     *     true
     */
    const degrees = get_degrees(graph);
    for (const key in degrees) {
        const sum = degrees[key].reduce((a, b) => a + b, 0);
        if (sum % 2 != 0) return false;
    }
    return true;
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

// BA03F: Eulerian cycle, by Hierholzer's algorithm
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

// BA03I: k universal circular string
function k_universal_circular_string(k) {
    /**
     * int -> str
     * returns a list of k_universal circular strings
     */
    const all_kmers = product_rec(['0', '1'], k);
    const graph = debruijn_from_kmers(all_kmers);
    if (!(check_eulerian_cycle(graph))) return "Not Eulerian!"
    const cycle = eulerian_hierholzer(graph);
    let answer = cycle[0];
    for (const elem of cycle.slice(1)) {
        answer += elem[elem.length - 1];
    }
    answer = answer.slice(0, answer.length - (k - 1));
    return answer;
}

// main
function main() {
    const fs = require('fs');
    // const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03i.txt').toString().split("\n");

    const startTime = performance.now()
    const cycle = k_universal_circular_string(4);
    console.log(cycle);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()
