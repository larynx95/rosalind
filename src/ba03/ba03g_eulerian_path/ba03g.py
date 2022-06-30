"""
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node) is called an Eulerian path.

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
      -> PREV: Eulerian Cycle (BA3F)
        -> NEXT: k-universal string problem (BA3I)
      -> HERE: Eulerian Path (BA3G)

Info.
* writing down graph to figure out what this exercise is:
  {0:[2], 1:[3], 2:[1], 3:[0, 4], 6:[3,7], 7:[8], 8:[9], 9:[6]}
                4
                ↑
        9 → 6 → 3 → 0       Not circular!!!
        ↑   ↓   ↑   ↓
        8 ← 7   1 ← 2
* problems:
  - This graph is not circular graph.
    Can any node be a starting point (source) or an ending point (target)?
  - How can I find the end node?
    a. comparing keys and values in dictionary (De Bruijn graph)
       dict.keys:   0123  6789
       dict.values: 012334 789   --> 4 is the end!
    b. difference bw in-degree and out-degree
       In-degree must be larger than out-degree in the last node.
       0:0, 1:0, 2:0, 3:0, 4:1, 6:-1, 7:0, 8:0, 9:0  --> The node 4 is the last!

Plan 1.
* re-using "EulerianCycle" function
  - Adding an extra edge from w (end node) to v (start node)
    transforms the Eulerian path into an Eulerian cycle.
  - convert Eulerian path to Eulerian cycle, then re-convert cycle to path

* out-degree +1, in-degree -1
                4
                ↑
        9 → 6 → 3 → 0
        ↑   ↓   ↑   ↓
        8 ← 7   1 ← 2
  unique nodes: 0    1    2    3    4    6    7    8    9
  --------------------------------------------------------
    0:[2]       +1       -1
    1:[3]           +1        -1
    2:[1]           -1   +1
    3:[0,4]     -1            +2   -1
    6:[3,7]                   -1        +2   -1
    7:[8]                                    +1   -1
    8:[9]                                         +1   -1
    9:[6]                               -1             +1
  --------------------------------------------------------
                0    0    0    0   -1   +1    0    0    0   --> start:6 end:4
* add edge (end, start) to graph, get Eulerian cycle: 4 -> 6
* trim the last repetitive elem in cycle, then rotate cycle to get Eulerian path

═════════════════════════════════════════════════

References:
- Remove all the elements that occur in one list from another
  https://stackoverflow.com/questions/4211209/remove-all-the-elements-that-occur-in-one-list-from-another
  >>> set([1,2,6,8]) - set([2,3,5,8])
      set([1, 6])
  >>> l3 = [x for x in l1 if x not in l2]
- How to make a flat list out of a list of lists?
  https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-a-list-of-lists
  >>> [item for sublist in lls for item in sublist]
- Is there a way I can initialize dictionary values to 0 in python taking keys from a list? [duplicate]
  https://stackoverflow.com/questions/52410880/is-there-a-way-i-can-initialize-dictionary-values-to-0-in-python-taking-keys-fro
"""
#!/usr/bin/env python
import time


def read_de_bruijn_graph_as_dict(lines):
    """
    De Bruijn Graph as dictionary
    >>> read_de_bruijn_graph_as_dict(lines)
        {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
    """
    adjlist = {}
    for line in lines:
        start, temp = [line.strip() for line in line.split(' -> ')]
        ends = [int(elem) for elem in temp.split(',')]
        adjlist[int(start)] = ends
    return adjlist


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


def eulerian_path_to_cycle(graph):
    """
    {a:[a]} -> {a:[a]}
    """
    (start,end) = find_src_target_nodes(graph)
    graph[end] = [start]
    return graph


def eulerian_cycle(graph):
    """
    {a:[a]} -> [a]
    returns Eulerian cycle
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


def eulerian_path(graph):
    """
    {a:[a]} -> [a]
    returns Eulerian path
    >>> eulerian_path({0:[2],1:[3],2:[1],3:[0,4],6:[3,7],7:[8],8:[9],9:[6]})
        [6,7,8,9,6,3,0,2,1,3,4]
    """
    # Eulerian path to Eulerian cycle
    (start,end) = find_src_target_nodes(graph)
    graph[end] = [start]
    cycle = eulerian_cycle(graph)
    # rotate
    cycle = cycle[:-1]
    idx_start = cycle.index(start)
    return cycle[idx_start:] + cycle[:idx_start]


def main():
    f = open('/home/wsl/rosalind/data/ba03g.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    graph = read_de_bruijn_graph_as_dict(lines)
    f.close()

    start_time = time.time()
    path = eulerian_path(graph)
    print('->'.join(list(map(str,path))))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()