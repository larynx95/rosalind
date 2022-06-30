"""
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
"""

import time
import sys


def de_bruijn(kmers):
    """
    [str] -> [(str:[str])]
    returns DeBruijnk graph from k-mer patterns
    representation of De Bruijn graph as list of tuples: not good idea, use dictionary
    >>> de_bruijn(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])
        [('AGG',['GGG']),('CAG',['AGG','AGG']),('GAG',['AGG']),('GGA',['GAG']),('GGG',['GGA','GGG'])]
    """
    dict = {}
    for kmer in kmers:
        key = kmer[:-1]
        val = kmer[1:]
        if key not in dict:
            dict[key] = [val]
        else:
            dict[key].append(val)
    lst = []
    for k, v in dict.items():
        lst.append((k,sorted(v)))  # list of tuples
    return sorted(lst)


def de_bruijn_dictionary(kmers):
    """
    [a] -> {a:[a]}
    returns DeBruijnk graph from k-mer patterns (BA3E)
    representation of De Bruijn graph as list of tuples: not good idea, use dictionary
    >>> de_bruijn(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])
        {'GAG':['AGG'],'CAG':['AGG','AGG'],'GGG':['GGG','GGA'],'AGG':['GGG'],'GGA':['GAG']}
    """
    result = {}
    for kmer in kmers:
        key = kmer[:-1]
        val = kmer[1:]
        if key not in result:
            result[key] = [val]
        else:
            result[key].append(val)
    return result


def print_de_bruijn(kmers):
    for elem in de_bruijn(kmers):
        print(elem[0], end=" -> ")
        print(*elem[1], sep =",")


def main():
    f = open('/home/wsl/rosalind/data/ba03e.txt', 'r')
    kmers = [kmer.strip() for kmer in f.readlines()]
    f.close()

    start_time = time.time()
    print_de_bruijn(kmers)
    print("--- %s seconds ---" % (time.time() - start_time))

    #start_time = time.time()
    #original_stdout = sys.stdout
    #with open('/home/wsl/rosalind/data/ba03e_output.txt', 'w') as f:
    #    sys.stdout = f
    #    print_de_bruijn(kmers)
    #    sys.stdout = original_stdout
    #print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()