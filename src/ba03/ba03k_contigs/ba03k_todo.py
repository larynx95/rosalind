"""
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

"""

import time


def allpaths(graph, start, sink, path=[]):
    """
    ({a:[a]},a,a,[a]) -> [[a]]
    returns all path from 'start' node to 'sink' node
    """
    path = path + [start]
    if start not in graph:
        return [path]
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = allpaths(graph, node, sink, path)
            for newpath in newpaths:
                if newpath[-1] == sink:
                    paths.append(newpath)
    return paths


def main():
    f = open('/home/wsl/rosalind/data/ba03k.txt', 'r')
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()