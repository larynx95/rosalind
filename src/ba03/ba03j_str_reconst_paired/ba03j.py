"""
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
    4 2          prefix     suffix         012345678901234567   PathGraph[i]
    -----------------------------------------------------------------------
    GTGG|GTGA -> (GTG|GTG)  (TGG|TGA)      GTGG__GTGA           0
    TGGT|TGAG -> (TGG|TGA)  (GGT|GAG)       TGGT__TGAG          1
    GGTC|GAGA -> (GGT|GAG)  (GTC|AGA)        GGTC__GAGA         2<-- gap was filled here
    GTCG|AGAT -> (GTC|AGA)  (TCG|GAT)         GTCG__AGAT        3<-- overlapping starts here
    TCGT|GATG -> (TCG|GAT)  (CGT|ATG)          TCGT__GATG       4
    CGTG|ATGT -> (CGT|ATG)  (GTG|TGT)           CGTG__ATGT      5
    GTGA|TGTT -> (GTG|TGT)  (TGA|GTT)            GTGA__TGTT     6
    TGAG|GTTG -> (TGA|GTT)  (GAG|TTG)             TGAG__GTTG    7
    GAGA|TTGA -> (GAG|TTG)  (AGA|TGA)              GAGA__TTGA   8
                                           GTGGTCGTGAGA
                                                 GTGAGATATTGA
    index:  012345678901
                  012345678901
    string1 GTGGTCGTGAGA
    string2       GTGAGATATTGA
            GTGGTCGTGAGATATTGA   == string1 + string2[k+d:]  TODO: Is this always right?

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
"""
#!/usr/bin/env python
import time


def paired_composition(text, k, d):
    """
    ([a],int,int) -> {([a],[a])}
    implementation of 'PairedComposition_(k,d)(Text)'
    returns a collection of 'PairedReads'
    >>> paired_composition('TAATGCCATG',3,1)
        [('TAA','GCC'),('AAT','CCA'),('ATG','CAT'),('TGC','ATG')]
    >>> paired_composition([0,1,2,3,4,5,6,7,8,9],3,1)
        [([0,1,2],[4,5,6]),([1,2,3],[5,6,7]),([2,3,4],[6,7,8]),([3,4,5],[7,8,9])]
    """
    result = []
    for i in range(len(text) - (2*k+d) + 1):
        pattern1 = text[i:i+k]
        pattern2 = text[i+k+d:i+2*k+d]
        result.append((pattern1, pattern2))
    return result


def composition_graph(text, k, d):
    """
    (str,int,int) -> [((str,str),(str,str))]
    implementation of 'CompositionGraph_(k,d)(Text)'
    returns a collection of isolated edges (list of tuples of two tuples)
    This is not useful! I can't reconstruct a string with this!
    >>> composition_graph('TAATGCCATG',3,1)
        [(('TA','GC'),('AA','CC')),(('AA','CC'),('AT','CA')),(('AT','CA'),('TG','AT')),(('TG','AT'),('GC','TG'))]
    >>> composition_graph('AAAAAAAA',3,1)  --> This is why I used a list of tuples of two tuples
        [(('AA','AA'),('AA','AA')),(('AA','AA'),('AA','AA'))]
    """
    result = []
    for i in range(len(text) - (2*k+d) + 1):
        pattern1 = text[i:i+k]
        pattern2 = text[i+k+d:i+2*k+d]
        prefix = (pattern1[:-d], pattern2[:-d])
        suffix = (pattern1[d: ], pattern2[d: ])
        result.append((prefix,suffix))
    return result


def de_bruijn(text, k, d):
    """
    (a,int,int) -> {(a,a):[(a,a)]}
    implementation of 'DeBruijn_(k,d)(Text)'
    """
    graph = {}
    for i in range(len(text) - (2*k+d) + 1):
        pattern1 = text[i:i+k]
        pattern2 = text[i+k+d:i+2*k+d]
        prefix = (pattern1[:-d], pattern2[:-d])
        suffix = (pattern1[d: ], pattern2[d: ])
        if not prefix in graph:
            graph[prefix] = [suffix]
        else:
            graph[prefix].append(suffix)
    return graph


def de_bruijn_from_paired_reads(reads):
    """
    [(a,a)] -> {(a,a):[(a,a)]}
    returns 'DeBruijnGraph' from a list of PairedReads
    """
    graph = {}
    for read in reads:
        pattern1 = read[0]
        pattern2 = read[1]
        prefix = (pattern1[:-1], pattern2[:-1])
        suffix = (pattern1[1: ], pattern2[1: ])
        if not prefix in graph:
            graph[prefix] = [suffix]
        else:
            graph[prefix].append(suffix)
    return graph


def eulerian_cycle(graph):
    """
    {b:[b]} -> [b]    <-- b :: (a,a)
    returns Eulerian cycle (BA3F)
    """
    cycle = [min(graph.keys())]
    while len(graph) > 0:
        if cycle[0] == cycle[-1]:
            while not cycle[0] in graph:
                cycle.pop(0)
                cycle.append(cycle[0])
        source = cycle[-1]
        cycle.append(graph[source].pop())
        if len(graph[source]) == 0: del graph[source]
    return cycle


def find_src_target_nodes(graph):
    """
    {a:[a]} -< (a,a)
    returns a tuple of nodes, (start, end)
    >>> find_src_target_nodes({0:[2],1:[3],2:[1],3:[0,4],6:[3,7],7:[8],8:[9],9:[6]})
        (6, 4)
    """
    # get a list of unique nodes
    keys = graph.keys()
    ls_values = graph.values()
    values = []
    for value in ls_values:
        values.extend(value)
    unique_nodes = set(keys).union(set(values))
    # get a list score: n(out-degrees) - n(in-degrees)
    scores = dict.fromkeys(unique_nodes,0)
    for key,values in graph.items():
        for value in values:
            scores[key] += 1
            scores[value] -= 1
    # find start, end nodes
    for key,value in scores.items():
        if value > 0:
            start = key
        if value < 0:
            end = key
    return (start, end)


def eulerian_path(graph):
    """
    {b:[b]} -> [b]    <-- b :: (a,a)
    returns Eulerian path (BA3G)
    """
    # Eulerian path to Eulerian cycle
    (start,end) = find_src_target_nodes(graph)
    graph[end] = [start]
    cycle = eulerian_cycle(graph)
    # rotate
    cycle = cycle[:-1]
    idx_start = cycle.index(start)
    return cycle[idx_start:] + cycle[:idx_start]


def path_graph(text, k, d):
    """
    (a,int,int) -> [(a,a)]    <-- a :: str
    implementation of 'PathGraph_(k,d)(Text)'
    >>> path_graph('TAATGCCATGGGATGTT',3,1)
        [('TA','GC'),('AA','CC'),('AT','CA'),('TG','AT'),('GC','TG'),('CC','GG'),
        ('CA','GG'),('AT','GA'),('TG','AT'),('GG','TG'),('GG','GT'),('GA','TT')]
    """
    dbgraph = de_bruijn(text, k, d)
    epath = eulerian_path(dbgraph)
    return epath


def path_graph_from_paired_reads(reads):
    """
    [(a,a)] -> [(a,a)]    <-- a :: str
    implementation of 'PathGraph_(k,d)(Text)'
    """
    dbgraph = de_bruijn_from_paired_reads(reads)
    epath = eulerian_path(dbgraph)
    return epath


def str_reconst_wrong(path, k, d):
    """
    [(a,a)] -> a
    returns a string from 'pathGraph_(k,d)(Text)'
    """
    string1 = path[0][0]
    string2 = path[0][1]
    for tup in path[1:]:
        string1 += tup[0][-1]
        string2 += tup[1][-1]
    return string1 + string2[k+d:]


def str_reconst(path):
    """
    [(a,a)] -> a
    returns a reconstructed string from the result of 'PathGraph_(k,d)(TGext)'
    """
    string1 = path[0][0]
    string2 = path[0][1]
    for tup in path[1:]:
        string1 += tup[0][-1]
        string2 += tup[1][-1]
    idx = 0
    for i in range(len(string1)):
        if not string2.find(string1[i:]) == -1:
            idx = i
            break
    return string1[:idx] + string2


def main():
    f = open('/home/wsl/rosalind/data/ba03j.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, d = [int(num.strip()) for num in lines[0].split()]
    paired_reads = []
    for line in lines[1:]:
        prefix, suffix = [rd.strip() for rd in line.split('|')]
        paired_reads.append((prefix, suffix))
    f.close()

    start_time = time.time()
    path = path_graph_from_paired_reads(paired_reads)
    print(str_reconst(path))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()