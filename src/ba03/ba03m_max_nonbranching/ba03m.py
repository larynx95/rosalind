"""
Rosalind: BA3M
Generate All Maximal Non-Branching Paths in a Graph

A node v in a directed graph Graph is called a 1-in-1-out node
if its indegree and outdegree are both equal to 1, i.e., in(v) = out(v) = 1.
We can rephrase the definition of a "maximal non-branching path"
from the main text as a path whose internal nodes are 1-in-1-out nodes
and whose initial and final nodes are not 1-in-1-out nodes.
Also, note that the definition from the main text does not handle the special case
when Graph has a connected component that is an isolated cycle,
in which all nodes are 1-in-1-out nodes.

The MaximalNonBranchingPaths pseudocode below generates all non-branching paths in a graph.
It iterates through all nodes of the graph that are not 1-in-1-out nodes
and generates all non-branching paths starting at each such node.
In a final step, MaximalNonBranchingPaths finds all isolated cycles in the graph.

╔═══════════════════════════════════════════════════════════════════════════════╗
║  MaximalNonBranchingPaths(Graph)                                              ║
║    Paths <- empty list                                                        ║
║    for each node v in Graph                                                   ║
║      if v is not a 1-in-1-out node                                            ║
║        if out(v) > 0                                                          ║
║          for each outgoing edge (v, w) from v                                 ║
║            NonBranchingPath <- the path consisting of the single edge (v, w)  ║
║              while w is a 1-in-1-out node                                     ║
║                extend NonBranchingPath by the outgoing edge (w, u) from w     ║
║                  w <- u                                                       ║
║              add NonBranchingPath to the set Paths                            ║
║    for each isolated cycle Cycle in Graph                                     ║
║      add Cycle to Paths                                                       ║
║    return Paths                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝

Maximal Non-Branching Path Problem
Find all maximal non-branching paths in a graph.

Given: The adjacency list of a graph whose nodes are integers.

Return: The collection of all maximal non-branching paths in the graph.

Sample Dataset
1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6

Sample Output
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7

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
    * NEXT: Read-Pair (BA3J)
                                         "
              not every Eulerian path in the paired de Bruijn graph
         constructed from a (k, d)-mer composition spells out a solution of
                 the String Reconstruction from Read-Pairs Problem.
                                         "
      -> PREV: Construct a String Spelled by a Gapped Genome Path (BA3L)
      -> HERE: Generate All Maximal Non-Branching Paths in a Graph (BA3M)

Plan 1.
  ╔═══════════════════════════════════════════════════════════════════════════════╗
  ║  MaximalNonBranchingPaths(Graph)                                              ║
  ║    Paths <- empty list                                                        ║
  ║    for each node v in Graph                                                   ║
  ║      if v is not a 1-in-1-out node                                            ║
  ║        if out(v) > 0                                                          ║
  ║          for each outgoing edge (v, w) from v                                 ║
  ║            NonBranchingPath <- the path consisting of the single edge (v, w)  ║
  ║              while w is a 1-in-1-out node                                     ║
  ║                extend NonBranchingPath by the outgoing edge (w, u) from w     ║
  ║                  w <- u                                                       ║
  ║              add NonBranchingPath to the set Paths                            ║
  ║    for each isolated cycle Cycle in Graph                                     ║
  ║      add Cycle to Paths                                                       ║
  ║    return Paths                                                               ║
  ╚═══════════════════════════════════════════════════════════════════════════════╝

  * sample dataset
    [1] -> 2 -> [3] -> [4]      6 -> 7 -> 6
                 ↓
                [5]
    v           v       v    <-- not 1-in-1-out nodes: 1,3,4,5
    1:[2] 2:[3] 3:[4,5] 6:[7] 7:[6]

  * non-branching paths
    1 -> 2 -> 3 ┌ 3 -> 4
                └ 3 -> 5
    6 -> 7 -> 6

  * degrees
    {1:[2], 2:[3], 3:[4,5], 6:[7], 7:[6]}            <-- graph
    {1:[0,1], 2:[0,1], 3:[0,2], 6:[0,1], 7:[0,1]}    <-- out-degrees: {source:[-, len(graph[source])]}
    {1:[0,1], 2:[0,1], 3:[1,2], 4:[1,0], 5:[1,0], 6:[1,1], 7:[1,1]}  <-- in-degrees

Plan 2.

"""
#!/usr/bin/env python
import time


def read_debruijn_graph(lines):
    """
    (str) -> {int:[int]}
    helper function for constructing Eulerian graph
    dictionary representation of De Bruijn Graph
    >>> read_debruijn_graph(lines)
        {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
    """
    adjlist = {}
    for line in lines:
        start, temp = [line.strip() for line in line.split(' -> ')]
        ends = [int(elem) for elem in temp.split(',')]
        adjlist[int(start)] = ends
    return adjlist


def get_degrees(graph):
    """
    {a:[a]} -> {a:[int,int]}
    returns a dictionary {node: [in-degrees, out-degrees]}
    >>> get_degrees({1:[2],2:[3],3:[4,5],6:[7],7:[6]})
        {1:[0,1],2:[1,1],3:[1,2],4:[1,0],5:[1,0],6:[1,1],7:[1,1]}
    """
    result = {}
    for source, targets in graph.items():
        if source not in result:
            result[source] = [0, len(targets)]
        else:
            result[source][1] += len(targets)
        for target in targets:
            if target not in result:
                result[target] = [1,0]
            else:
                result[target][0] += 1
    return result


def remove_edge(graph, source, target):
    """
    ({a:[a]},a,a) -> {a:[a]}
    returns a modified graph in which an edge was removed
    """
    graph[source].remove(target)
    if not graph[source]:
        del graph[source]
    return graph


def max_nonbrancing_path(graph):
    """
    {a:[a]} -> [[a]]
    returns a list of lists of non-branching paths
    using textbook algorithm
    """
    paths = []
    degrees = get_degrees(graph)
    for source in degrees.keys():
        if degrees[source] != [1, 1]:
            if degrees[source][1] > 0:
                while source in graph:
                    target = graph[source][0]
                    non_branching_path = [source, target]
                    graph = remove_edge(graph, source, target)
                    while degrees[target] == [1, 1]:
                        next_target = graph[target][0]
                        non_branching_path.append(next_target)
                        graph = remove_edge(graph, target, next_target)
                        target = next_target
                    paths.append(non_branching_path)
    while graph:
        start = list(graph.keys())[0]
        cur = graph[start][0]
        graph = remove_edge(graph, start, cur)
        cycle = [start, cur]
        while cur != start:
            target = graph[cur][0]
            cycle.append(target)
            graph = remove_edge(graph, cur, target)
            cur = target
        paths.append(cycle)
    return paths


def main():
    f = open('/home/wsl/rosalind/data/ba03m.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    graph = read_debruijn_graph(lines)
    f.close()

    start_time = time.time()
    lls =max_nonbrancing_path(graph)
    for ls in lls:
        print(' -> '.join(list(map(str, ls))))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()