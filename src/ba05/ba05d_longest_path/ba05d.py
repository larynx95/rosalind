"""
Rosalind: BA5D
Find the Longest Path in a DAG

Longest Path in a DAG Problem
Find a longest path between two nodes in an edge-weighted DAG.

Given:
An integer representing the "source node" of a graph,
followed by an integer representing the "sink node" of the graph,
followed by an edge-weighted graph.
The graph is represented by a modified adjacency list in which the notation "0->1:7"
indicates that an edge connects node 0 to node 1 with weight 7.

Return:
The length of a longest path in the graph, followed by a longest path.
(If multiple longest paths exist, you may return any one.)

Sample Dataset
0
4
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3

Sample Output
9
0->2->3->4

═════════════════════════════════════════════════

Info.
  * representations for graph with weighted-edges (TODO: Which is the best?)
    a. list of triple: [(source, target, weight)]
       [(0, 1, 7), (0, 2, 4), (2, 3, 2), (1, 4, 1), (3, 4, 3)]
    b. dictionary: {(source, target): weight}
       {(0,1):7, (0,2):4, (2,3):2, (1,4):1, (3,4):3}
    c. dictionary: {source: [(target, weight)]}
       {0:[(1,7),(2,4)], 2:[(3,2), 1:[(4,2)], 3:[(4,3)]}
    d. dictionary: {source: ([edge],[weight])}
       {0:([1,2],[7,4]), 2:([3],[2]), 1:([4],[1]), 3:([4],[3])}

  * Dynamic programming in an arbitrary DAG
    (Important: provided that graph is topologically ordered)

    S(b) =            max             { S(a) + weight of edge from a to b }
          all predecessors a of nod b

  * In rectangular graph with diagonal edges

                 ┌ s(0,1) + weight of edge ↓  connecting (0,1) to (1,1)
    S(1,1) = max ├ s(1,0) + weight of edge →  connecting (1,0) to (1,1)
                 └ s(0,0) + weight of edge ↘ connecting (0,0) to (1,1)

                 ┌ s(i-1,j)   + weight of edge ↓  connecting (i-1,j) to (i,j)
    S(i,j) = max ├ s(i,j-1)   + weight of edge →  connecting (i,j-1) to (i,j)
                 └ s(i-1,j-1) + weight of edge ↘ connecting (i-1,j-1) to (i,j)

Plan 1. (Failed)
  * Pseudocode:
    ╔════════════════════════════════════════════════════════════════════════════════════╗
    ║ LONGESTPATH(Graph, source, sink)                                                   ║
    ║     for each node b in Graph                                                       ║
    ║         s_b <- -infinite                                                           ║
    ║     s_source <- 0                                                                  ║
    ║     topologically order Graph                                                      ║
    ║     for each node b in Graph (following the topological order)                     ║
    ║         s_b <- max<all predecessors a of node b>{s_a + weight of edge from a to b} ║
    ║     return s_sink                                                                  ║
    ╚════════════════════════════════════════════════════════════════════════════════════╝

  * steps: (This steps is not appropriate to solve this exercise.)
    (1) find topological order
    (2) find "start" and "sink" nodes
      - find nodes with no in-degree edges (start nodes)
      - find nodes with no out-degree edges (sink nodes)
    (3) construct dynamic table
      - assign 0 to all start nodes in dynamic table
      - from start nodes, construct dynamic table with a list which is topologically ordered in step (1)
      - at the end of this step, each node has its max length
      - I can get the the length of longest path of a node by indexing dynamic table. "table[node]"
    (4) find path
      - from the sink nodes, backtrack path to the start nodes
      - at the end of this step, I can get one or more longest paths

  * sample dataset
    - diagram
      edges : [(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9]

                    3
                ┌────────┐          4      2      3
        4    2  |     1  ↓      [0] -> [2] -> [3] -> [4]
      0 -> 2 -> 3   1 -> 4       0      4      6      9
      │            ↑↑               7      1                  9      1
      └────────────┘│9          [0] -> [1] -> [4]         [5] -> [1] - > [4]
           7        5            0      7      8           0      9      10

      S(4) = max(S(3)+3, S(1)+1)
                 ┬───    ┬───
                 │       └ max(S(0)+7, S(5)+9)
                 │             ┬───    ┬───
                 │             0       0
                 S(2)+2
                 ┬───
                 └ S(0)+4

    - constructing dynamic table
      edges                                              dynamic table (as dictionary)
      --------------------------------------------------------------------------------
                            topological ordering  -->     0   5   1   2   3   4
      {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}  {0:0,5:0}
        b       b                               b         *   *
      {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}  {0:0,5:0,1:9}
        a b                                     a b           --- *
      {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}  {0:0,5:0,1:9,2:4}
                a b                                       ---         *
      {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}  {0:0,5:0,1:9,2:4,3:6}
                     ^  a b                                           --- *
      {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}  {0:0,5:0,1:9,2:4,3:6,4:10}  <-- final table
                               a1 b    a2 b          ^            ---         *

    - get longest path
      {0:0,5:0,1:9,2:4,3:6,4:10}   {(0,1):7,(0,2):4,(2,3):2,(1,4):1,(3,4):3,(5,1):9}
      dynamic table, path          in-degree edges  max calculation
      --------------------------------------------------------------------------------
      table[4] == 10               max ┌ (1,4):1   table[1] + edges[(1,4)] = 9+1 = 10 *
                                       └ (3,4):3   table[3] + edges[(3,4)] = 6+3 = 9
      table[1] == 9                max ┌ (0,1):7   table[0] + edges[(0,1)] = 0+7 = 7
                                       └ (5,1):9   table[5] + edges[(5,1)] = 0+9 = 9 *
      table[5] == 0                end

Plan 2.
  * Plan 1 was failed.
  * In this exercise ...
    - The dataset given in this exercise is too simple to properly grasp the problem.
    - Two nodes were already given to me. (not always 'start' and 'sink' nodes)
    - Path may not exist, only one may exist, or more than one may exist.

  * steps:
    (1) find topological order --> a list of nodes
      - There can be many topological ordering.
    (2) find paths from the given "start" node  --> [(node, length)]
      - topo:   9     3     ...
        path: [(9:0),(3,4), ... ]
      - paths
        i. each path starts from the given start node and ends at the given sink node
        ii. path with longest length is the answer

Plan 3.
  * get a sublist from a list of topological ordered nodes (start:9, sink:19)
      I don't need all the nodes in topologically ordered list.
      ['2','5','0','1','3','11','7','10','4','9','6','12','8','14','15','13','16','18','19','17']
                                         => ['9','6','12','8','14','15','13','16','18','19']
      ['2','5','11','0','1','6','3','4','12','9','7','10','8','14','15','17','13','18','16','19']
                                         => ['9','7','10','8','14','15','17','13','18','16','19']
      ['2','0','4','1','5','11','3','6','8','14','7','13','15','10','17','9','18','12','16','19']
                                                                     => ['9','18','12','16','19']
      But ... which one? => Anyone will do cuz every list is already topologically ordered.(TODO: Is it right?)

  * get paths
    - example
        ┌────────────────────────────────────────────┐┌───────────┐
        │┌──↓┌────↓                                  ↓│           ↓
        9    7    10    8    14    15    17    13    18    16    19
        │                          ↑│           │          ↑│    ↑↑
        │                          ││           └──────────┘└────┘│
        └──────────────────────────┘└─────────────────────────────┘
    - candidates
      path 1: 9 -> 18 -> 19
      path 2: 9 -> 7 -> 10         --> should be discarded
      path 3: 9 -> 15 -> 19
      path 4: 13 -> 16 -> 19       --> not required
      others:                      --> not required
    - How to get all candidate-paths
      Can I get all paths by modifying BA3M solution (get all non-branching paths)?
    - choosing correct path
      Does the path always start at the start node and end at the sink node? => No.
      How can I choose correct one?
        a. start at node "start, end at node "sink"
        b. path with max length

Plan 4.
  * steps:
    - get all paths from 'start' node to 'end' node <--- DFS, BFS
    - find a path with the longest length

═════════════════════════════════════════════════

References:
- How to initialize a dict with keys from a list and empty value in Python?
  https://stackoverflow.com/questions/2241891/how-to-initialize-a-dict-with-keys-from-a-list-and-empty-value-in-python
- Dynamic programming vs Backtracking
  https://www.javatpoint.com/dynamic-programming-vs-backtracking
- Dynamic programming: top-down vs. bottom-up
  - What is the difference between bottom-up and top-down?
    https://stackoverflow.com/questions/6164629/what-is-the-difference-between-bottom-up-and-top-down
- BA3M (Is this pseudocode useful?)
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
- Find all paths in DAG from starting node **
  https://stackoverflow.com/questions/54441582/find-all-paths-in-dag-from-starting-node
"""

import time
import random


#################################################
# (failed) topological ordering => longest path TODO: Fix this
#################################################


def get_start_nodes(graph):
    """
    {int:[int]} -> {int}
    returns a set of start nodes
    >>> get_start_nodes({0:[1,2], 2:[3], 3:[4], 4:[5]})
        {0}
    """
    fsts = set()
    snds = set()
    for key, value in graph.items():
        fsts.add(key)
        snds |= set(value)
    return fsts - snds


def get_target_nodes(node_start, graph):
    """
    (int, {int:[int]}) -> {int}
    returns a set of target nodes
    >>> get_target_nodes(0,{0:[1,2],2:[3],3:[4],4:[5]})
        {1,2}
    """
    nodes_end = set()
    for key, value in graph.items():
        if key == node_start:
            nodes_end |= set(value)
    return nodes_end


def has_indegree(node, graph):
    """
    (int,{int:[int]}) -> bool
    returns True if a node has at least one in-degree edge
    >>> has_indegree(0,{0:[1,2],2:[3],3:[4],4:[5]})
        False
    >>> has_indegree(2,{0:[1,2],2:[3],3:[4],4:[5]})
        True
    """
    targets = set()
    for key, value in graph.items():
        targets |= set(value)
        if node in targets:
            return True
    return False


def get_indegree_edges(node, edges):
    """
    (a,{(a,a):int}) -> [(a,a)]
    returns a list of in-degree edges of a specific node
    >>> get_indegree_edges('1',{('0','1'):7,('0','2'):4,('2','3'):2,('1','4'):1,('3','4'):3,('5','1'):9})
        [('0','1'),('5','1')]
    >>> get_indegree_edges('1',{('0','2'):7})
        []
    """
    indegrees = []
    for edge, weight in edges.items():
        if edge[1] == node:
            indegrees.append((edge[0], edge[1]))
    return indegrees


def topological_ordering(graph):
    """
    {int:[int]} -> [int]
    returns topologically ordered graph (BA5N), algorithm in textbook
    >>> topological_ordering({0:[1,2],2:[3],3:[4],4:[5]})
        [0,2,3,4,1,5]
    """
    result = []
    candidates = get_start_nodes(graph)
    while candidates:
        a_node = random.choice(list(candidates))
        candidates.remove(a_node)
        result.append(a_node)
        for b_node in get_target_nodes(a_node, graph):
            graph[a_node].remove(b_node)
            if len(graph[a_node]) == 0:
                del graph[a_node]
            if not has_indegree(b_node, graph):
                candidates.add(b_node)
    if graph:
        return "The input graph is not a DAG."
    else:
        return result


def longest_path_wrong(edges):
    """
    {(a,a):int} -> (int, [a])
    returns the longest path
    incorrect - In this fx, source and sink nodes were not specified.
    TODO: Don't rush. Read the exercises thoroughly before typing.
    """
    # edges to graph without weights
    graph = {}
    for tup in edges.keys():
        source, target = [tup[0], tup[1]]
        if source not in graph:
            graph[source] = [target]
        else:
            graph[source].append(target)
    # construct dynamic table as dictionary
    table = {}
    starts = get_start_nodes(graph)
    for start in starts:
        table[start] = 0
    ordered = topological_ordering(graph)
    for node in ordered:
        if node in table:
            continue
        max_value = 0
        for edge in get_indegree_edges(node, edges):
            val = table[edge[0]] + edges[edge]
            if val > max_value:
                max_value = val
        table[node] = max_value
    # get the length of longest path
    max_length = 0
    for node, length in table.items():
        if length > max_length:
            max_length = length
    # find path
    sinks = []
    for node, weight in table.items():
        if weight == max_length:
            sinks.append(node)
    paths = []
    for sink in sinks:
        path = [sink]
        target = sink
        while table[target] != 0:
            tups = get_indegree_edges(target, edges)
            max_val = 0
            source = ''
            for tup in tups:
                val = table[tup[0]] + edges[(tup[0], target)]
                if val > max_val:
                    max_val = val
                    source = tup[0]
            path.append(source)
            target = source
        paths.append(path[::-1])
    return max_length, paths


#################################################
# all paths => longest path
#################################################


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


def get_length(edges, path):
    """
    ({(a,a):int},{a:[a]}) -> int
    returns length of a path
    >>> get_length({('0','1'):7,('0','2'):4,('2','3'):2,('1','4'):1,('3','4'):3,('5','1'):9}, ['0','1','4'])
        8
    """
    length = 0
    for i in range(len(path)-1):
        length += edges[(path[i], path[i+1])]
    return length


def longest_path(edges, start, end):
    """
    {(a,a):int} -> (int, [a])
    returns the longest length, path
    """
    # edges to graph without weights
    graph = {}
    for tup in edges.keys():
        source, target = [tup[0], tup[1]]
        if source not in graph:
            graph[source] = [target]
        else:
            graph[source].append(target)
    # get all path from 'start' to 'end' node
    paths = allpaths(graph, start, end)
    # find a path with the longest length
    max_length = 0
    longest_path = []
    for path in paths:
        leng = get_length(edges, path)
        if  leng > max_length:
            max_length = leng
            longest_path = path
    return max_length, longest_path


def main():
    with open('/home/wsl/rosalind/data/ba05d.txt', 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        start = lines[0]   # start node - A node doesn't need to be a integer.
        end = lines[1]    # sink node
        edges = {}         # edges :: {(a,a):int}
        graph = {}         # graph :: {a:[a]}
        for line in lines[2:]:
            fragments = [frag for frag in line.split(':')]
            weight = int(fragments[1].strip())
            source, target = [elem.strip() for elem in fragments[0].split('->')]
            # edges
            edge = (source,target)
            edges[edge] = weight
            # graph
            if source not in graph:
                graph[source] = [target]
            else:
                graph[source].append(target)

    # sample data for practice
    sample_edges = {('0','1'): 7,('0','2'): 4,('2','3'): 2,('1','4'): 1,('3','4'): 3,('5','1'): 9}
    sample_graph = {'0':['1','2'], '1':['4'], '2':['3'], '3':['4'], '5':['1']}

    # execution
    start_time = time.time()
    print(allpaths(sample_graph, '0','4'))
    answers = longest_path(edges, start, end)
    print(answers[0])
    print('->'.join(answers[1]))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
