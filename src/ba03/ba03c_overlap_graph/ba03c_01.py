'''
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

Sample Output
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG

-------------------------------------------------

solution by Elmar Hinz
'''

import time
from collections import defaultdict

def formatGraph(graph):
    lines = []
    for source in sorted(graph):
        for target in sorted(graph[source]):
            lines.append("%s -> %s" % (source, target))
    return ("\n").join(lines)

def OverlapGraph(nodes):
   # 1. group target nodes by prefix
    targets = defaultdict(list)
    for node in nodes:
        prefix = node[:-1]
        targets[prefix].append(node)
    # 2. map target nodes to postfix
    graph = {}
    for node in nodes:
        postfix = node[1:]
        if postfix in targets:
            graph[node] = targets[postfix]
    return graph

def main():
    f = open('/home/wsl/rosalind/data/ba03c.txt', 'r')
    patterns = [line.strip() for line in f.readlines()]
    f.close()

    start_time = time.time()
    print(formatGraph(OverlapGraph(patterns)))
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()