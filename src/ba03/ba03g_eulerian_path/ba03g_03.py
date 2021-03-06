'''
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node is called an Eulerian path.

Eulerian Path Problem
Find an Eulerian path in a graph.

Given: A directed graph that contains an Eulerian path,
       where the graph is given in the form of an adjacency list.

Return: An Eulerian path in this graph.

Sample Dataset
0 -> 2
1 -> 3
2 -> 1
3 -> 0,4
6 -> 3,7
7 -> 8
8 -> 9
9 -> 6

Sample Output
6->7->8->9->6->3->0->2->1->3->4

-------------------------------------------------

https://github.com/dohlee/daily-rosalind/blob/master/scripts/bioinformatics-textbook-track/BA3G.py
'''

##################################################
# Find an Eulerian Path in a Graph
#
# http://rosalind.info/problems/BA3G/
#
# Given: A directed graph that contains an Eulerian
#  path, where the graph is given in the form of
#  an adjacency list.
#
# Return: An Eulerian path in this graph.
#
# AUTHOR : dohlee
##################################################

# Your imports here
from BA3F import parse_adjacency_list, find_eulerian_cycle
from collections import Counter

# Your codes here
def add_imaginary_edge(graph):
    """Add imaginary edge (start, end), where start is the node
    with surplus incoming edges, and end is the node with surplus
    outgoing edges.
    """
    outgoingEdgeCounts, incomingEdgeCounts = Counter(), Counter()
    for u in graph:
        outgoingEdgeCounts[u] += len(graph[u])
        for v in graph[u]:
            incomingEdgeCounts[v] += 1

    start = list((incomingEdgeCounts - outgoingEdgeCounts).keys())[0]
    end = list((outgoingEdgeCounts - incomingEdgeCounts).keys())[0]

    # Add imaginary edge.
    graph[start].append(end)
    return start, end

def find_eulerian_path(graph):
    """Return an eulerian path of the graph."""
    start, end = add_imaginary_edge(graph)
    cycle = find_eulerian_cycle(graph, start=end)[:-1]
    for i in range(len(cycle)):
        if cycle[i] == start and cycle[(i+1) % len(cycle)] == end:
            path = cycle[i+1:] + cycle[:i+1]

    return path

if __name__ == '__main__':
    # Load the data.
    with open('/home/wsl/rosalind/data/ba03g.txt') as inFile:
        graph = parse_adjacency_list(inFile)

    # Print output
    with open('/home/wsl/rosalind/data/ba03g_output.txt', 'w') as outFile:
        print('->'.join(map(str, find_eulerian_path(graph))), file=outFile)