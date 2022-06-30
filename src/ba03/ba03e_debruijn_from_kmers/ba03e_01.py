'''
Rosalind: BA3E
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

-------------------------------------------------

by apense
https://rosalind.info/users/apense/
'''

import time
import sys

def de_bruijn_list(dna_list):
    result = {}
    for dna in dna_list:
        try:
            result[dna[:-1]].append(dna[1:])
        except KeyError:
            result[dna[:-1]] = [dna[1:]]
    return result

def main():
    f = open('/home/wsl/rosalind/data/ba03e.txt', 'r')
    kmers = [kmer.strip() for kmer in f.readlines()]
    f.close()

    start_time = time.time()
    m = de_bruijn_list(kmers)
    original_stdout = sys.stdout
    with open('/home/wsl/rosalind/data/ba03e_output.txt', 'w') as f:
        sys.stdout = f
        for key in m:
            f.write(key + " -> " + ','.join([q for q in m[key]]) + "\n")
        sys.stdout = original_stdout
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()