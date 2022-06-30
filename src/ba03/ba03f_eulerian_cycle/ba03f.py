"""
Rosalind: BA3F (difficulty: 4/5)
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

╔══════════════════════════════════════════════════════════════════════╗
║ EULERIANCYCLE(Graph)                                                 ║
║     form a cycle Cycle by randomly walking in Graph                  ║
║         (don't visit the same edge twice!)                           ║
║     while there are unexplored edges in Graph                        ║
║         select a node newStart in Cycle with still unexplored edges  ║
║         form Cycle' by traversing Cycle (starting at newStart)       ║
║           and then randomly walking                                  ║
║         Cycle <- Cycle'                                              ║
║     return Cycle                                                     ║
╚══════════════════════════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║
║         v <- a node with indegree larger than 1 in G        ║
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║
║                 if NewGraph is connected                    ║
║                     add NewGraph to ALlGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝

Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An Eulerian directed graph, in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3
1 -> 0
2 -> 1,6
3 -> 2
4 -> 2
5 -> 4
6 -> 5,8
7 -> 9
8 -> 7
9 -> 6

Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

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
      -> PREV: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> HERE: Eulerian Cycle (BA3F)
      -> NEXT: Eulerian Path (BA3G)

Info.
- analyzing sample dataset:

      4 ← 5
      ↓   ↑            In cycle, any node (vertex) can be a start node.
  1 ←[2]→[6]→ 8        Every node has equal number of in- and out-degree.
  ↓   ↑   ↑   ↓        If a node has multiple out-degrees, I should select correct one among them.
  0 → 3   9 ← 7

- What is the "Eulerian Cycle" and what conditions make a graph Eulerian?
  ; every vertex has even degree
  ; every vertex has equal in degree and out degree
  ; all of its vertices with nonzero degree belong to a single connected component (??)

- balanced graph:
  A vertex j in V is said to be balanced if its indegree and the outdegree are equal.
  Now the graph is said to be balanced if all the vertices are balanced.

- De Bruijn Graph --> Eulerian path
  Overlap Graph   --> Hamiltonian path

- representations of De Bruijn Graph
  - a list of lists       :  [[source, target...]]
  - a collection of tuples:  [(source,target)] or {(source,target)}
  - dictionary            :  {source:[target]}

- walk simulation
  Eulerian Graph                                                     next cycle
  ------------------------------------------------------------------------------
  {1:[0],0:[3],3:[2],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  1    []
  {1:[-],0:[3],3:[2],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  0    [1]
  {1:[-],0:[-],3:[2],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  3    [1,0]
  {1:[-],0:[-],3:[-],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  2    [1,0,3]
  {1:[-],0:[-],3:[-],2:[-,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  1    [1,0,3,2]         dead end, wrong selection
  {1:[-],0:[-],3:[-],2:[1,-],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  6    [1,0,3,2]         return to the past
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  5    [1,0,3,2,6]       wrong selection
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,8],5:[-],4:[2],8:[7],7:[9],9:[6]}  4    [1,0,3,2,6,5]
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,8],5:[-],4:[2],8:[7],7:[9],9:[6]}  2    [1,0,3,2,6,5,4]
  {1:[-],0:[-],3:[-],2:[-,-],6:[-,8],5:[-],4:[2],8:[7],7:[9],9:[6]}  1    [1,0,3,2,6,5,4,2] dead end
  {1:[-],0:[-],3:[-],2:[1,-],6:[5,-],5:[4],4:[2],8:[7],7:[9],9:[6]}  8    [1,0,3,2,6]       return to the past
  {1:[-],0:[-],3:[-],2:[1,-],6:[5,-],5:[4],4:[2],8:[-],7:[9],9:[6]}  7    [1,0,3,2,6,8]
  {1:[-],0:[-],3:[-],2:[1,-],6:[5,-],5:[4],4:[2],8:[-],7:[-],9:[6]}  9    [1,0,3,2,6,8,7]
  {1:[-],0:[-],3:[-],2:[1,-],6:[5,-],5:[4],4:[2],8:[-],7:[-],9:[-]}  6    [1,0,3,2,6,8,7,9]
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,-],5:[4],4:[2],8:[-],7:[-],9:[-]}  5    [1,0,3,2,6,8,7,9]
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,-],5:[-],4:[2],8:[-],7:[-],9:[-]}  4    [1,0,3,2,6,8,7,9,5]
  {1:[-],0:[-],3:[-],2:[1,-],6:[-,-],5:[-],4:[-],8:[-],7:[-],9:[-]}  2    [1,0,3,2,6,8,7,9,5,4]
  {1:[-],0:[-],3:[-],2:[-,-],6:[-,-],5:[-],4:[-],8:[-],7:[-],9:[-]}  1    [1,0,3,2,6,8,7,9,5,4,2]
  {1:[-],0:[-],3:[-],2:[-,-],6:[-,-],5:[-],4:[-],8:[-],7:[-],9:[-]}  -    [1,0,3,2,6,8,7,9,5,4,2,1]

- Problems:
  a. How do I know that I made the wrong choice at a node with multiple out-degrees?
  b. How can I go back to the point when I made the wrong choice?

Plan 1.
  - find any generalized rule: There're multiple small cycles (or circles)!

          4 ← 5        1->0->3->[2]  stop!
          ↓   ↑        1->0->3→2     →1→x dead end
      1 ←[2]→[6]→ 8                  →6   →5→4→2→1→x dead end
      ↓   ↑   ↑   ↓                       →8→7→9→6→5→4→2→1  success
      0 → 3   9 ← 7    In the left graph, node '2' has only two out-degrees.
                      But what happens if node '2' has a lot of out-degrees?

    if a node has a out-degree, then just walk through the Eulerian graph until ...
      a. there's another node with multiple edges
      b. there's no unexplored edge

  - divide Eulerian graph into small circles, then merge those small circles into a large cycle
    - selecting start node: node with minimum value
    - merging circles: TODO: How can I merge all small circles into one large circle?
      This is not easy at all.
      There is one common node between two small circles.
      It is not eacy to find out which circles have common node.(!!!)

    A. case 1
          4 <- 5                          4 <- 5
          ↓    ↑                          ↓    ↑
      1 <-[2]->[6]-> 8        1 <-[2]    [2]->[6]    [6]-> 8
      ↓    ↑    ↑    ↓        ↓   ↑                  ↑    ↓
      0 -> 3    9 <- 7        0 -> 3                  9 <- 7
                              circle1    circle2     circle3
      circles           unexplored edges
      -------------------------------------------
      0->3->2->1->0->x  {2:[6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  circle1
      2->6->5->4->2->x  {6:[8],8:[7],7:[9],9:[6]}                      circle2
      6->8->7->9->6->x  {}                                             circle3

    B. case 2
          4 <- 5                         4 <- 5
          ↓    ↑                         ↓    ↑
      1 <-[2]->[6]-> 8        1 <-[2]    [2]->[6]-> 8
      ↓    ↑    ↑    ↓        ↓   ↑           ↑    ↓
      0 -> 3    9 <- 7        0 -> 3           9 <- 7
                              circle1     circle2
      circles                       unexplored edges
      -------------------------------------------
      0->3->2->1->0->x              {2:[6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}  circle1
      2->6->8->7->9->6->5->4->2->x  {}                                             circle2

Plan 2.
  - start from small circle, rotate it, then move to the next small circle
    Key point: walk --> a circle --> rotate --> walk --> ...

           4 <- 5
           ↓    ↑             counterclockwise (or cloclwise) rotation
      1 <-[2]->[6]-> 8        0 <- 1     1 <- 2
      ↓    ↑    ↑    ↓        ↓   ↑  ==  ↓   ↑
      0 -> 3    9 <- 7        3 -> 2     0 -> 3

    circle                                 closed
    -----------------------------------------------------------
    0->3->2->1->0                          Yes
    3->2->1->0->3                          r
    2->1->0->3->2                          r  node '2' has remaining edge(s)
    2->1->0->3->2->6->5->4->2              Yes
    1->0->3->2->6->5->4->2->1              r
    0->3->2->6->5->4->2->1->0              r
    3->2->6->5->4->2->1->0->3              r
    2->6->5->4->2->1->0->3->2              r
    6->5->4->2->1->0->3->2->6              r  node '6' has remaining edge(s)
    6->5->4->2->1->0->3->2->6->8->7->9->6  complete

Plan 3.
- mysterious algorithm, clojure (by hadaarjan)
- Sample algorithm by hadaarjan (Brilliant!)
  He (or She) focused on the recursive aspect of this probem.
  Structure:
  * Outside the recursion
    - a list for storing path
    - a starting node
  * base case (condition for ending recursion)
    - if selected key has no value
  * recursive case
  * visualization
    Example A.
                 ↙ ↖
      0 → 1 →[2]→ 3 → 4
        ↖ ↙

    def eulerianCycle_hadaarjan(gdict):
        start = min(gdict.keys())
        tour = []
        def cycle(n):
            while len(gdict[n]) > 0:
                next = gdict[n].pop()
                cycle(next)
            tour.append(n)
        cycle(start)
        tour = tour[::-1]
        return tour

    graph                              call stack                               tour
    ---------------------------------------------------------------------------------
    {0:[1],1:[2],2:[3,0],3:[4],4:[2]}  f(0)                                      []
    {0:[-],1:[2],2:[3,0],3:[4],4:[2]}  f(0)->f(1)                                []
    {0:[-],1:[-],2:[3,0],3:[4],4:[2]}  f(0)->f(1)->f(2)                          []
    {0:[-],1:[-],2:[3,-],3:[4],4:[2]}  f(0)->f(1)->f(2)->f(0)                    []
    {0:[-],1:[-],2:[3,-],3:[4],4:[2]}  f(0)->f(1)->f(2)                          [0]
    {0:[-],1:[-],2:[-,-],3:[4],4:[2]}  f(0)->f(1)->f(2)->f(3)                    [0]
    {0:[-],1:[-],2:[-,-],3:[-],4:[2]}  f(0)->f(1)->f(2)->f(3)->f(4)              [0]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)->f(1)->f(2)->f(3)->f(4)->f(2)        [0]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)->f(1)->f(2)->f(3)->f(4)              [02]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)->f(1)->f(2)->f(3)                    [024]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)->f(1)->f(2)                          [0243]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)->f(1)                                [02432]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  f(0)                                      [024321]
    {0:[-],1:[-],2:[-,-],3:[-],4:[-]}  the end                                   [0243210]

    >>> graph = {0:[1],1:[2],2:[3,0],3:[4],4:[2]}
    >>> eulerianCycle_hadaarjan(graph)
        n tour               graph
    -------------------------------------------------------------------
    A:  0 []                 {0:[1],1:[2],2:[3,0],3:[4],4:[2]}
    B:  0 1                  {0:[ ],1:[2],2:[3,0],3:[4],4:[2]}
    A:  1 []                 {0:[ ],1:[2],2:[3,0],3:[4],4:[2]}
    B:  1 2                  {0:[ ],1:[ ],2:[3,0],3:[4],4:[2]}
    A:  2 []                 {0:[ ],1:[ ],2:[3,0],3:[4],4:[2]}
    B:  2 0                  {0:[ ],1:[ ],2:[3  ],3:[4],4:[2]}
    A:  0 []                 {0:[ ],1:[ ],2:[3  ],3:[4],4:[2]}
    C:  0 [0]                {0:[ ],1:[ ],2:[3  ],3:[4],4:[2]}
    B:  2 3                  {0:[ ],1:[ ],2:[   ],3:[4],4:[2]}
    A:  3 [0]                {0:[ ],1:[ ],2:[   ],3:[4],4:[2]}
    B:  3 4                  {0:[ ],1:[ ],2:[   ],3:[ ],4:[2]}
    A:  4 [0]                {0:[ ],1:[ ],2:[   ],3:[ ],4:[2]}
    B:  4 2                  {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    A:  2 [0]                {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  2 [0,2]              {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  4 [0,2,4]            {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  3 [0,2,4,3]          {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  2 [0,2,4,3,2]        {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  1 [0,2,4,3,2,1]      {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    C:  0 [0,2,4,3,2,1,0]    {0:[ ],1:[ ],2:[   ],3:[ ],4:[ ]}
    -------------------------------------------------------------------
    result: [0,1,2,3,4,2,0]

Plan 4.
- no rotation, no merging, just stop and remember state

           4 <- 5
           ↓    ↑
      1 <-[2]->[6]-> 8
      ↓    ↑    ↑    ↓
      0 -> 3    9 <- 7

  circle
  -------------------------------------------------------
  0->3->2  stop! remember current state! next node: 1,6
           ->1->0->x  dead end, discard
           ->6        stop! remember current state! next node: 5,8
                      ->5->4->2->1->0->x  dead end, discard
                      ->8->7->9->6->5->4->2->1->0  complete! ... smell of recursion ...
  TODO: implement this function with recursion (interesting!)

- simulation:

  circle  state
          walk with copied graph    remaining unexplored edges
  ---------------------------------------------------------------------------------------
  032     {1:[0],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}
          032-10                     {      2:[  6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}
          032-6                      {1:[0],2:[1  ],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}
  0326    {1:[0],2:[1  ],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]}
          0326-54210                 {              6:[  8],            8:[7],7:[9],9:[6]}
          0326-879654210             {}

- stop condition:
  a. no out-degree with no remaining unexplored edge --> success
  b. no out-degree with remaining unexplored edge    --> failed, dead end
  c. node with just one out-degree (edge)            --> walk
  d. node with two or more out-degree (edges)        --> iterate wakling for all edges

Plan 5.
- Non-Binary Tree data structure, depth first search (DFS) algorithm, consider OOP
  TODO: try OOP later
  [(0,3),(1,0),(2,1),(2,6),(3,2),(4,2),(5,4),(6,5),(6,8),(7,9),(8,7),(9,6)]
  {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  (0,[3])          (0,3)        0:[3]
  (1,[0])          (1,0)        1:[0]
  (2,[1,6])        (2,1)        2:[1,6]
                   (2,6)
  (3,[2])          (3,2)        3:[2]
  (4,[2])     or   (4,2)  or    4:[2]
  (5,[4])          (5,4)        5:[4]
  (6,[5,8])        (6,5)        6:[5,8]
                   (6,8)
  (7,[9])          (7,9)        7:[9]
  (8,[7])          (8,7)        8:[7]
  (9,[6])          (9,6)        9:[6]

      03
      |
      32
     /  \
   21    26
   |    /  \
  10  65    68
      |     |
      54    87
      |     |
      42    79
      |     |
      21    96
      |     |
      10    65
            |
            54
            |
            42
            |
            21
            |
            10

═════════════════════════════════════════════════

References:
  * Basics
    - Eulerian path
      https://en.wikipedia.org/wiki/Eulerian_path
    - signed graph
      https://en.wikipedia.org/wiki/Signed_graph
    - Does Eulerian cycle in digraph really need strongly connected component?
      https://math.stackexchange.com/questions/2296452/does-eulerian-cycle-in-digraph-really-need-strongly-connected-component
    - Euler’s Theorem.
      If a pseudograph G has an Eulerian circuit, then G is connected and the degree of every vertex is even.
    - Hierholzer’s Theorem.
      If a pseudograph G is connected and the degree of every vertex of G is even, then G has an Eulerian circuit.
      If every vertex in a pseudograph G has positive even degree, then any given vertex of G lies on some circuit of G.
    - Veblen’s Theorem.
      A pseudograph G has a decomposition into cycles if and only if every vertex of G has even degree.
    - Listing’s Theorem.
      If G is a connected pseudograph with precisely sh vertices of odd degree, h 6= 0,
      then there exists h trails in G such that each edge of G is exactly on of these trails.
      Furthermore, fewer than h trails with this property cannot be found.
  * Hierholzer's algorithm
    - Eulerian Path/Circuit algorithm (Hierholzer's algorithm) | Graph Theory
      https://www.youtube.com/watch?v=8MpoO2zA2l4
    - Existence of Eulerian Paths and Circuits | Graph Theory
      https://www.youtube.com/watch?v=xR4sGgwtR2I
    - Example of Hierholzer's algorithm
      https://www.youtube.com/watch?v=qZrfK2iE4UA
  * Fleury's Algorithm
    - Eulerian Path and Circuit Algorithm - How does it work?
      https://cs.stackexchange.com/questions/145496/eulerian-path-and-circuit-algorithm-how-does-it-work
    - Euler Part 3: Fleury's Algorithm for Finding an Euler Circuit in Graph with Vertices of Even Degree
      https://www.youtube.com/watch?v=F4BM6fnLl04
  * Veblen's theorem
      https://en.wikipedia.org/wiki/Veblen%27s_theorem
  * bridge finding algorithms
    - cut (graph theory)
      https://en.wikipedia.org/wiki/Cut_(graph_theory)
    - Intro to Menger's Theorem | Graph Theory, Connectivity
      https://www.youtube.com/watch?v=u1BiIP0zw-c
    - Eunice’s algorithm - Determine if edge is a bridge in a graph
      https://math.stackexchange.com/questions/3965493/determine-if-edge-is-a-bridge-in-a-graph
    - Find Bridges in a graph using Tarjans Algorithm | Cut Edge
      https://www.youtube.com/watch?v=Rhxs4k6DyMM
    - chain decompositions, Dilworth's theorem
      https://en.wikipedia.org/wiki/Dilworth%27s_theorem
      What are Graph Decompositions? | Graph Decomposition, Graph Theory
      https://www.youtube.com/watch?v=BiLjdZYj_RI
    - Tarjans Strongly Connected Components algorithm | Graph Theory
      https://www.youtube.com/watch?v=TyWtx7q2D7Y
    - Kosaraju's algorithm
      https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm
      Strongly Connected Components Kosaraju's Algorithm Graph Algorithm
      Kosaraju algorithm introduction
      https://www.youtube.com/watch?v=Jb1XlDsr46o
  * Deep copy of a dict in python
    https://stackoverflow.com/questions/5105517/deep-copy-of-a-dict-in-python
  * Finding the index of an item in a list
    https://stackoverflow.com/questions/176918/finding-the-index-of-an-item-in-a-list
  * Testing equality of three values
    https://stackoverflow.com/questions/13805939/testing-equality-of-three-values
  * Usage of the "==" operator for three objects
    https://stackoverflow.com/questions/13792604/usage-of-the-operator-for-three-objects/13792615#13792615
"""
#!/usr/bin/env python
import time


###########################################################
# Hierholzer's algorithm
###########################################################


def Hierholzer(graph):
    """
    {int:[int]}
    returns a Eulerian path or cycle by Hierholzer algorithm
    """
    pass


###########################################################
# Fleury's Algorithm
###########################################################


def Fleury(graph):
    """
    {int:[int]}
    retunrs a Eulerain path or cycle by Fleury's Algorithm
    """
    pass


###########################################################
# Plan 1. merging small circles (Veblen's theorem??) TODO: Failed1 Fix this!
###########################################################


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


def walk_dirty(graph):
    """
    {a:[a]} -> [a]
    >>> walk_dirty({1:[0],0:[3],3:[2],2:[1,6],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [0,3,2,6,5,4,2,1,0]
    >>> walk_dirty({1:[0],0:[3],3:[2],2:[6,1],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [0,3,2,1,0]
    """
    cur = min(graph.keys())
    path = [cur]
    while True:
        if cur not in graph:           # if a node does not in graph
            break
        next_nodes = graph[cur]        # get values (list) of key 'cur'
        next_node = next_nodes.pop()   # pop last value
        if next_nodes == []:           # if values is [] after popping, delete key
            del graph[cur]
        path.append(next_node)         # append a node to path
        cur = next_node                # update cur node
    return path


def remove_edge(adj_list, from_node, to_node):
    """
    ({a:[a]},a,a) -> {a:[a]}
    """
    adj_list[from_node].remove(to_node)
    if not adj_list[from_node]:
        del adj_list[from_node]
    return adj_list


def walk(graph):
    """
    {a:[a]} -> [a]
    return a circle, clean version
    >>> walk({1:[0],0:[3],3:[2],2:[1,6],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [0,3,2,1,0]
    >>> walk({1:[0],0:[3],3:[2],2:[6,1],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [0,3,2,6,8,7,9,6,5,4,2,1,0]
    """
    cur = min(graph.keys())
    path = [cur]
    while cur in graph:
        # appending a node to a circle
        next_nodes = graph[cur]
        next_node = next_nodes[0]
        path.append(next_node)
        # remove current edge
        remove_edge(graph,cur,next_node)  # <-- simplified external function
        # update current node
        cur = next_node
    return path


def get_circles(graph):
    """
    {a:[a]} -> [a]
    returns a list of small circles from Eulerian graph
    TODO: Well... Can I merge small circles? How?
    >>> get_circles({1:[0],0:[3],3:[2],2:[1,6],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [[0,3,2,6,5,4,2,1,0],[6,8,7,9,6]]
    >>> get_circles({1:[0],0:[3],3:[2],2:[6,1],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [[0,3,2,1,0],[2,6,5,4,2],[6,8,7,9,6]]
    >>> get_circles({1:[0],0:[3],3:[2],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [[0,3,2,6,8,7,9,6,5,4,2,1,0]]
    """
    circles = []
    while graph:
        circle = walk(graph)
        circles.append(circle)
    return circles


def merge_two_circles_wrong(circ1, circ2):
    """
    ([a],[a]) -> [a]
    return merged circle
    TODO: incorrect! This function only works when each node has at most two edges.
    >>> merge_two_circles_wrong([0,3,2,1,0],[2,6,5,4,2])
        [2,1,0,3,2,6,5,4,2]
    """
    common_node = list(set(circ1).intersection(set(circ2)))[0]  # TODO: only ONE?
    idx1 = circ1.index(common_node)
    idx2 = circ2.index(common_node)
    return circ1[:-1][idx1:] + circ1[:-1][:idx1] + circ2[:-1][idx2:] + circ2[:-1][:idx2] + [circ2[idx2]]


def merge_circles_wrong(circles):
    """
    [[a]] -> [a]
    merge small circles into one large circle
    TODO: incorrect! This function only works when each node has at most two edges.
    >>> merge_circles_wrong([[0,3,2,1,0],[6,5,4,2,6],[6,8,7,9,6]])
        [6,5,4,2,1,0,3,2,6,8,7,9,6]
    """
    result = circles[0]
    for circle in circles[1:]:
        result = merge_two_circles_wrong(result, circle)
    return result


def eulerian_cycle_wrong(graph):
    """
    {a:[a]} -> [a]
    return Eulerian cycle by merging small circles
    TODO: incorrect! This function only works when each node has at most two edges.
    >>> eulerian_cycle_wrong({1:[0],0:[3],3:[2],2:[1,6],6:[8,5],5:[4],4:[2],8:[7],7:[9],9:[6]})
        [6,5,4,2,1,0,3,2,6,8,7,9,6]
    """
    circles = get_circles(graph)
    return merge_circles_wrong(circles)


###########################################################
# Plan 2. Rotation and walk
###########################################################


def eulerian_cycle(graph):
    """
    solution by Elmar Hinz
    very clever!
    """
    cycle = [min(graph.keys())]                       # init cycle with any node
    while len(graph) > 0:                             # run until all edges are moved from graph to cycle
        if cycle[0] == cycle[-1]:                     # whenever cycle closes, rotate to a node with remaining targets
            while not cycle[0] in graph:
                cycle.pop(0)
                cycle.append(cycle[0])
        source = cycle[-1]                            # the last node of cycle is the new source
        cycle.append(graph[source].pop())             # move one target at a time from graph to the end of cycle
        if len(graph[source]) == 0: del graph[source] # clean up empty dict entries of graph
    return cycle


###########################################################
# Plan 3. Mysterious algorithm
###########################################################


def eulerianCycle_hadaarjan(gdict):
    """
    solution by hadaarjan
    """
    start = min(gdict.keys())
    tour = []
    def cycle(n):                   #  ┐ nested function
        while len(gdict[n]) > 0:    #  │
            next = gdict[n].pop()   #  │
            cycle(next)             #  │
        tour.append(n)              #  ┘
    cycle(start)
    tour = tour[::-1]
    return tour


###########################################################
# Plan 4. TODO: Recursive
###########################################################


###########################################################
# Plan 5. TODO: Non-Binary Tree, OOP
###########################################################


###########################################################
# TODO: All Eluerian Cycles
###########################################################


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


def is_simple_graph(graph):
    """
    {a:[a]} -> bool
    returns True if a graph is simple graph (in-degrees == out-degrees == 1)
    """
    degrees = get_degrees(graph)
    for node in degrees:
        if not node[0] == node[1] == 1:
            return False
    return True


def all_nonsimple_graphs(ls_graph):
    """
    [{a:[a]}] -> bool
    returns a list of all non-simple graphs
    """
    graphs = []
    for graph in ls_graph:
        if not is_simple_graph(graph):
            graphs.append(graph)
    return graphs


def node_with_indegrees(graph):
    """
    {a:[a]} -> a
    returns a node with in-degree more than 1
    """
    degrees = get_degrees(graph)
    for key, value in degrees:
        if value[0] > 1:
            return key


def get_incomming_edges(graph, vertex):
    """
    ({a:[a]},a) -> [(a,a)]
    returns all incomming edges
    >>> get_incomming_edges({1:[0],0:[3],3:[2],2:[1,6],6:[5,8],5:[4],4:[2],8:[7],7:[9],9:[6]},2)
        [(3,2),(4,2)]
    """
    edges = []
    for source, targets in graph.items():
        for target in targets:
            if vertex == target:
                edges.append((source, target))
    return edges


def get_outgoing_edges(graph, vertex):
    """
    ({a:[a]},a) -> [(a,a)]
    returns all outgoing edges
    """
    return graph[vertex]


def bypass(graph, u, v, w):
    """
    ({a:[a]},a,a,a) -> {a:[a]}
    returns (u,v,w)-bypass graph - fig.ba3f.3.41
    TODO: What exactly is the "bypass"?
    """
    pass


def all_eulerian_cycles(graph):
    """
    {a:[a]} -> {[a]}
    """
    allgraphs = [graph]
    while True:
        nonsimple_graphs = all_nonsimple_graphs(allgraphs)
        if not nonsimple_graphs:
            break
        for nsgraph in nonsimple_graphs:
            vertex = node_with_indegrees(nsgraph)
            outgoing_edges  = get_outgoing_edges(nsgraph, vertex)
            for in_edge in get_incomming_edges(nsgraph, vertex):
                for out_edge in get_outgoing_edges(nsgraph, vertex):
                    pass # TODO: complete this function

"""
╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║
║         v <- a node with indegree larger than 1 in G        ║
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║
║                 if NewGraph is connected                    ║
║                     add NewGraph to ALlGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝
"""


def main():
    f = open('/home/wsl/rosalind/data/ba03f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    graph = read_debruijn_graph(lines)

    #start_time = time.time()
    #cycle = eulerian_cycle_wrong(graph)
    #for elem in cycle:
    #  print(elem, end=" ")
    #print()
    #print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    cycle = eulerian_cycle(graph)
    print(cycle)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
