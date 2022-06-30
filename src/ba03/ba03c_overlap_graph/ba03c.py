"""
Rosalind: BA3C
Construct the Overlap Graph of a Collection of k-mers

In this chapter, we use the terms prefix and suffix to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.

Given an arbitrary collection of k-mers Patterns, we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern'). The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

Overlap Graph Problem
Construct the overlap graph of a collection of k-mers.

Given: A collection Patterns of k-mers.

Return: The overlap graph Overlap(Patterns), in the form of an adjacency list.

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
- How to modify list entries during for loop?
  https://stackoverflow.com/questions/4081217/how-to-modify-list-entries-during-for-loop
- Sorting a set of values
  https://stackoverflow.com/questions/17457793/sorting-a-set-of-values
- Index all *except* one item in python
  https://stackoverflow.com/questions/19286657/index-all-except-one-item-in-python
- How to deep copy a list?
  https://stackoverflow.com/questions/17873384/how-to-deep-copy-a-list
"""

import time


def overlap_graph_adjlist(patterns):
    """ Overlap(Patterns) """
    starts = patterns[:]
    result = []
    for pat in starts:
        for target in patterns:
            if pat[1:] == target[:-1]:
                result.append((pat, target))
    return result


def print_adjlist(adjlist):
    sorted_adjlist = sorted(adjlist)
    for elem in sorted_adjlist:
        print(elem[0] + " -> " + elem[1])


def main():
    f = open('/home/wsl/rosalind/data/ba03c.txt', 'r')
    patterns = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    adjlist = overlap_graph_adjlist(patterns)
    print_adjlist(adjlist)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()