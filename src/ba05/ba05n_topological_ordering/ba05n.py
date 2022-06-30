"""
Rosalind: BA5N
Find a Topological Ordering of a DAG

Topological Ordering Problem
Find a topological ordering of a directed acyclic graph.

Given: The adjacency list of a graph (with nodes represented by integers).

Return: A topological ordering of this graph.

Sample Dataset (Be careful! There may be multiple target nodes: "1 -> 2,3,4,5")
1 -> 2
2 -> 3
4 -> 2
5 -> 3

Sample Output
1, 4, 5, 2, 3

═════════════════════════════════════════════════

    [ Where am I? ]

    * learn dynamic programming (DP) in "coin change problem" (BA5A)
      ↓
    * apply dynamic programming to "Manhattan tourist problem" (BA5B)
      ↓
    * PREV: apply solution of "Manhattan tourist problem" to "longest common subsequence (LCS) problem" (BA5C)
      ↓
    * NEXT: get generalized version of algorithm for any Directed Acyclic Graph (DAG) (BA5D)
        <- HERE: Topological ordering (BA5N)


Plan 1.

- How can I represent adjacent list?
  a. list of tuples : [(0,1), (1,5), ... ]
  b. dictionary     : {0:[1,2,3], 2:[6,7,8], ... }  <-- use this!

- What is the "topological ordering"?
  an ordering of nodes (a1, . . . , ak) in a DAG is called a topological ordering
  if every edge (ai, aj) of the DAG connects a node
  with a smaller index to a node with a larger index, i.e., i < j.

- topological ordering in practice, some cases (algorithm in textbook):
  ╔════════════════════════════════════════════════════════════════════╗
  ║ TOPOLOGICALORDERING(Graph)                                         ║
  ║     List <- empty list                                             ║
  ║     Candidates <- set of all nodes in Graph with no incoming edges ║
  ║     while Candidates is non-empty                                  ║
  ║         select an arbitrary node a from Candidates                 ║
  ║         add a to the end of List and remove it from Candidates     ║
  ║         for each outgoing edge from a to another node b            ║
  ║             remove edge (a, b) from Graph                          ║
  ║             if b has no incoming edges                             ║
  ║                 add b to Candidates                                ║
  ║     if Graph has edges that have not been removed                  ║
  ║         return "the input graph is not a DAG"                      ║
  ║     else                                                           ║
  ║         return List                                                ║
  ╚════════════════════════════════════════════════════════════════════╝

  practice 1:
    ┌────────┐                            ┌──┐┌────────┐
    │        ↓                            │  ↓|        ↓
    0 -> 1   2 -> 3 -> 4             0 -> 2   3   1 -> 4
         │             ↑             │            ↑
         └─────────────┘             └────────────┘
    Graph                             a  (a,b)     Candidates      list
    -----------------------------------------------------------------
    [(0,1),(0,2),(2,3),(1,4),(3,4)]                    {0}         []
    [  -     -   (2,3),(1,4),(3,4)]   0  (0,1),(0,2)   {}          [0]
    [  -     -   (2,3),(1,4),(3,4)]                    {1,2}              <-- select an arbitrary node a from Candidates: choose 1, not 2
    [  -     -   (2,3)   -   (3,4)]   1  (1,4)         {2}         [0,1]
    [  -     -     -     -   (3,4)]   2  (2,3)         {3}         [0,1,2]
    [  -     -     -     -     -  ]   3  (3,4)         {4}         [0,1,2,3]
    [  -     -     -     -     -  ]   4                {}          [0,1,2,3,4]

    practice 2:
         ┌──┐┌────────┐
         │  ↓|        ↓
    0 -> 2   3   1 -> 4
    │            ↑
    └────────────┘
    Graph                             a  (a,b)     Candidates      list
    -----------------------------------------------------------------
    [(0,1),(0,2),(2,3),(1,4),(3,4)]                    {0}         []
    [  -     -   (2,3),(1,4),(3,4)]   0  (0,1),(0,2)   {}          [0]
    [  -     -   (2,3),(1,4),(3,4)]                    {1,2}              <-- choose 2, not 1
    [  -     -     -   (1,4),(3,4)]   2  (2,3)         {1,3}       [0,2]
    [  -     -     -   (1,4)   -  ]   3  (3,4)         {1}         [0,2,3]
    [  -     -     -      -    -  ]   1  (1,4)         {4}         [0,2,3,1]
    [  -     -     -      -    -  ]   4                {}          [0,2,3,1,4]

    What I found out in above two practices:
    - Weigths of edges are not required for topological ordering
    - There can be multiple topological ordering.
    - In DAG, there should be a source node and a sink node.

    practice 3:
    - multiple source nodes! --> change the plan! use dictionary as adjacent list, and modify code

          4   5
          ↓   ↓
      1 → 2 → 3

═════════════════════════════════════════════════

References:
- Topological sorting
  https://en.wikipedia.org/wiki/Topological_sorting#:~:text=In%20computer%20science%2C%20a%20topological,before%20v%20in%20the%20ordering.
  In computer science,
  a topological sort or topological ordering of a directed graph
  is a linear ordering of its vertices
  such that for every directed edge uv from vertex u to vertex v,
  u comes before v in the ordering.
- In a DAG, does there exist a sink for all vertices in the DAG?
  https://cs.stackexchange.com/questions/67125/in-a-dag-does-there-exist-a-sink-for-all-vertices-in-the-dag
- Can a DAG have multiple source nodes? (yes)
  https://www.iarcs.org.in/inoi/online-study-material/topics/dags.php
- How to retrieve an element from a set without removing it?
  https://stackoverflow.com/questions/59825/how-to-retrieve-an-element-from-a-set-without-removing-it
- Retrieve elements from a Python set
  https://www.techiedelight.com/retrieve-elements-from-set-python/
- random.choice from set? python
  https://stackoverflow.com/questions/15837729/random-choice-from-set-python
- How to join two sets in one line without using "|"
  https://stackoverflow.com/questions/17429123/how-to-join-two-sets-in-one-line-without-using
- Python: Checking if a 'Dictionary' is empty doesn't seem to work
  https://stackoverflow.com/questions/23177439/python-checking-if-a-dictionary-is-empty-doesnt-seem-to-work
- Delete an element from a dictionary
  https://stackoverflow.com/questions/5844672/delete-an-element-from-a-dictionary
- other pseudocodes:
╔════════════════════════════════════════════════════════════════════════════════════════╗
║ MERGE(List1, List2)                                                                    ║
║     SortedList <- empty list                                                           ║
║     while both List1 and List2 are non-empty                                           ║
║         if the smallest element in List1 is smaller than the smallest element in List2 ║
║             move the smallest element from List1 to the end of SortedList              ║
║         else                                                                           ║
║             move the smallest element from List2 to the end of SortedList              ║
║     move any remaining elements from either List1 or List2 to the end of SortedList    ║
║     return SortedList                                                                  ║
╚════════════════════════════════════════════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════════════════╗
║ MERGESORT(List)                                             ║
║     if List consists of a single element                    ║
║         return List                                         ║
║     FirstHalf <- first half of List                         ║
║     SecondHalf <- second half of List                       ║
║     SortedFirstHalf <- MERGESORT(FirstHalf)                 ║
║     SortedSecondHalf <- MERGESORT(SecondHalf)               ║
║     SortedList <- MERGE(SortedFirstHalf, SortedSecondHalf)  ║
║     return SortedList                                       ║
╚═════════════════════════════════════════════════════════════╝

"""
#!/usr/bin/env python3
import time
import random


#################################################
# list of tuples version (has limitations, not good)
#################################################


def find_start_nodes_from_adjlist(adjlist):
    """
    [(int,int)] -> {int}
    find source and sink nodes from adjacent list (list of tuples)
    source node: node wihtout in-degree
    sink node  : node without out-degree
    >>> find_start_nodes_from_adjlist([(0,1),(0,2),(2,3),(1,4),(3,4)])
        {0}
    >>> find_start_nodes_from_adjlist([(0,1),(0,2),(5,2),(2,3),(1,4),(3,4)])
        {0,5}
    """
    fsts = set()
    snds = set()
    for tup in adjlist:
        fsts.add(tup[0])
        snds.add(tup[1])
    return fsts - snds


def find_target_nodes_from_adjlist(node_start, adjlist):
    """
    (a,[(a,a)]) -> {a}
    returns end nodes from start node 'node_start'
    >>> find_target_nodes_from_adjlist(0,[(0,1),(0,2),(2,3),(1,4),(3,4)])
        {1, 2}
    """
    nodes_end = set()
    for tup in adjlist:
        if tup[0] == node_start:
            nodes_end.add(tup[1])
    return nodes_end


def has_indegree_in_adjlist(node, adjlist):
    """
    (a,[(a,a)]) -> bool
    returns True if a node has indegree
    >>> has_indegree_in_adjlist(0,[(2,3),(1,4),(3,4)])
        False
    >>> has_indegree_in_adjlist(4,[(2,3),(1,4),(3,4)])
        True
    """
    for tup in adjlist:
        if node == tup[1]:
            return True
    return False


def topological_ordering_from_adjlist(adjlist):
    """
    [(a,a)] -> [a]
    implemented as similar to the algorithm of the text book as possible
    >>> topological_ordering_from_adjlist([(0,1),(0,2),(2,3),(1,4),(3,4)])
        [0,2,1,3,4]    # <-- [0,1,2,3,4],[0,2,3,1,4],...
    """
    result = []
    candidates = find_start_nodes_from_adjlist(adjlist)
    while candidates:
        node_a = random.choice(list(candidates))
        candidates.remove(node_a)
        result.append(node_a)
        for node_b in find_target_nodes_from_adjlist(node_a, adjlist):
            adjlist.remove((node_a, node_b))
            if not has_indegree_in_adjlist(node_b, adjlist):
                candidates.add(node_b)
    if not len(adjlist) == 0:
        return "The input graph is not a DAG."
    else:
        return result


#################################################
# dictionary version (more general, and better)
#################################################


def find_start_nodes_from_adjdict(adjdict):
    """
    {int:[int]} -> {int}
    >>> find_start_nodes_from_adjdict({0:[1,2], 2:[3], 3:[4], 4:[5]})
        {0}
    """
    fsts = set()
    snds = set()
    for key, value in adjdict.items():
        fsts.add(key)
        snds |= set(value)
    return fsts - snds


def find_target_nodes_from_adjdict(node_start, adjdict):
    """
    (int, {int:[int]}) -> {int}
    >>> find_target_nodes_from_adjdict(0,{0:[1,2],2:[3],3:[4],4:[5]})
        {1,2}
    """
    nodes_end = set()
    for key, value in adjdict.items():
        if key == node_start:
            nodes_end |= set(value)
    return nodes_end


def has_indegree_in_adjdict(node, adjdict):
    """
    (int,{int:[int]}) -> bool
    >>> has_indegree_in_adjdict(0,{0:[1,2],2:[3],3:[4],4:[5]})
        False
    >>> has_indegree_in_adjdict(2,{0:[1,2],2:[3],3:[4],4:[5]})
        True
    """
    targets = set()
    for key, value in adjdict.items():
        targets |= set(value)
        if node in targets:
            return True
    return False


def topological_ordering_from_adjdict(adjdict):
    """
    {int:[int]} -> [int]
    >>> topological_ordering_from_adjdict({0:[1,2],2:[3],3:[4],4:[5]})
        [0, 2, 3, 4, 1, 5]
    """
    result = []
    candidates = find_start_nodes_from_adjdict(adjdict)
    while candidates:
        a_node = random.choice(list(candidates))
        candidates.remove(a_node)
        result.append(a_node)
        for b_node in find_target_nodes_from_adjdict(a_node, adjdict):
            adjdict[a_node].remove(b_node)
            if len(adjdict[a_node]) == 0:
                del adjdict[a_node]
            if not has_indegree_in_adjdict(b_node, adjdict):
                candidates.add(b_node)
    if adjdict:
        return "The input graph is not a DAG."
    else:
        return result


def main():
    try:
        with open('/home/wsl/rosalind/data/ba05n.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            adjdict = {}
            for line in lines:
                temp = [part.strip() for part in line.split('->')]
                start = int(temp[0])
                ends = [int(elem) for elem in temp[1].split(',')]
                adjdict[start] = ends
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    answer = topological_ordering_from_adjdict(adjdict)
    answer = list(map(str, answer))
    print(', '.join(answer))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
