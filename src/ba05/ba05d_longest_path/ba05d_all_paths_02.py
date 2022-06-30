"""
Topic: get all paths from a graph
  1. Representation of graph
    1) List vs. Set vs. Dictionary
      - Set (or dictionary) is better than list type because there's no priority or order in edges.
      - Dictionary has many advantages over Set type.
      - Set of edges :: {(a,a)}
      - Dictionary of edges :: {a:[a]} or :: {(a,a):int}
    2) Class vs. set (or dictionary)
      - Trie, Tree, ...
  2. various paths
    - get all paths from an arbitrary start node to an arbitrary sink node
    - get all paths from a start node            to an arbitrary node
    - get all paths from an arbitrary node       to a sink node
    - get all paths from a start node            to a sink node

═════════════════════════════════════════════════

Example 1:
  edges = [(6, 7), (6, 15), (15, 16), (16, 21), (15, 9), (9, 13), (13, 4), (4, 1), (1, 5)]
  graph = {6:[7,15], 15:[9,16], 16:[21], 9:[13], 13:[4], 4:[1], 1:[5]}

   6 -> 7
   ↓
  15 -> 9 -> 13 -> 4 -> 1 -> 5
   ↓
  16 -> 21

Output 2:

  6->7
  6->15->9->13->4->1->5
  6->15->16->21

═════════════════════════════════════════════════

Example 2:

  graph = {'a':['b'],'d':['e'],'k':['i'],'q':['r'],'c':['f'],'g':['h'],'j':['i'],'b':['c','e'],'e':['f','j'],'f':['g'],'i':['f','l'],'s':[],'t':[],'m':['n'],'n':['o'],'o':['p'],'p':['m'],'v':['w'],'w':['x'],'x':['n'] }

      a -> b -> c                p -> m       s      t
           ↓    ↓                ↑    ↓
      d -> e -> f -> g -> h      o <- n -> v
           ↓    ↑                     ↑    ↓
           j -> i <- k                x <- w
                ↓
                l                q -> r

Output 2:
  [['a', 'b', 'c', 'f', 'g', 'h'],
   ['a', 'b', 'e', 'f', 'g', 'h'],
   ['a', 'b', 'e', 'j', 'i', 'f', 'g', 'h'],
   ['a', 'b', 'e', 'j', 'i', 'l']]

═════════════════════════════════════════════════

References:
- Find all paths in DAG from starting node **
  https://stackoverflow.com/questions/54441582/find-all-paths-in-dag-from-starting-node
"""

import time


def paths1(graph, start, path=[]):
    path = path + [start]
    if start not in graph:
        return [path]
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = paths1(graph, node, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def paths2(graph, start, path=[]):
    path = path + [start]
    paths = []    # <-- `paths = [path]` is wrong
    if start not in graph:
        return [path]
    for node in graph[start]:
        newpaths = paths2(graph, node, path)
        for newpath in newpaths:
            paths.append(newpath)
    return paths


def paths3(graph, start, path=[]):
    path = path + [start]
    paths = [path]
    if start not in graph:
        return paths
    for node in graph[start]:
        newpaths = paths3(graph, node, path)
        for newpath in newpaths:
            paths.append(newpath)
    return paths


def paths4(graph, start):
    if start not in graph:
        return [[start]]
    paths = []
    for v in graph[start]:
        for v_path in paths4(graph, v):
            paths.append([start] + v_path)
    return paths


def main():
    # data
    edges = [(6, 7), (6, 15), (15, 16), (16, 21), (15, 9), (9, 13), (13, 4), (4, 1), (1, 5)]
    graph = {6:[7,15], 15:[9,16], 16:[21], 9:[13], 13:[4], 4:[1], 1:[5]}

    # execution
    start_time = time.time()
    print("paths 1: ", paths1(graph, 6))
    print("paths 2: ", paths2(graph, 6))
    print("paths 3: ", paths3(graph, 6))
    print("paths 4: ", paths4(graph, 6))
    graph = {'a':['b'],'d':['e'],'k':['i'],'q':['r'],'c':['f'],'g':['h'],'j':['i'],'b':['c','e'],'e':['f','j'],'f':['g'],'i':['f','l'],'s':[],'t':[],'m':['n'],'n':['o'],'o':['p'],'p':['m'],'v':['w'],'w':['x'],'x':['n'] }
    print(paths4(graph, 'e'))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
