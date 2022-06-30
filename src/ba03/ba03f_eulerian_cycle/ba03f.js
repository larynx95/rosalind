/*
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

Given: An "Eulerian directed graph", in the form of an adjacency list.

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

═════════════════════════════════════════════════

Info.
  * What is the "Eulerian Cycle" and what conditions make a graph Eulerian?
    ; every vertex has even degree
    ; every vertex has equal in degree and out degree
    ; all of its vertices with nonzero degree belong to a single connected component (??)

  * balanced graph:
    A vertex j in V is said to be balanced if its indegree and the outdegree are equal.
    Now the graph is said to be balanced if all the vertices are balanced.

  * De Bruijn Graph --> Eulerian path
    Overlap Graph   --> Hamiltonian path

  * representations of De Bruijn Graph
    - a list of lists       :  [[source, target...]]
    - a collection of tuples:  [(source,target)] or {(source,target)}
    - dictionary            :  {source:[target]}

  * walk simulation
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

  * Problems (motivation):
    a. How do I know that I made the wrong choice at a node with multiple out-degrees?
    b. How can I go back to the point when I made the wrong choice?

═════════════════════════════════════════════════ Algorithms

1. Hierholzer's algorithm
  * steps:
    ; just walk through the graph, and get a circle1
    ; if the circle1 is complete whole large cycle without remaining nodes in graph, that's Eulerian cycle
    ; if there're remaining unvisited node in graph yet, just do another walking, get another cycle2
      but in this time, start node should be a node in previous cycle, and it has one or more out-degrees
    ; patch cycle2 to cycle1, repeat previous steps
              4 < 5             4 < 5             4 * 5             4 * 5
              v   ^             v   ^             *   *             *   *
          1 < 2 > 6 > 8     1 *[2]> 6 > 8     1 * 2 * 6 > 8     1 * 2 *[6]* 8
          v   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
          0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
                            10321             10321             103265421
                                                 26542              68796
                                              103265421         1032687965421
  * another steps (very complex, not recommended):
    ; construct a I/O table of In/Out degrees
    ; find start node with Degree(out) - Degree(in) = 1
      if no such node, graph is a Eulerian cycle or it is not Eulerian
    ; focus only on Degree(out) in I/O table
    ; walk from the start node, delete edges, reduce the number of Degree(out) in the I/O table
    ; when the number of target node's Degree(out) in I/O table is zero, and any remaining node exists, move back to the previous node
    ; continue untill there's no edge in a graph
              4 < 5             4 < 5             4 < 5             4 < 5             4 * 5             4 < 5             4 * 5
              v   ^             v   ^             v   ^             v   ^             *   *             v   ^             *   *
          1 < 2 > 6 > 8     1 * 2 > 6 > 8     1 < 2 > 6 > 8     1 < 2 * 6 > 8     1 * 2 * 6 > 8     1 < 2 * 6 * 8     1 * 2 * 6 * 8
          v   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
          0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
    node: 0123456789        0123456789        0123456789        0123456789        0123456789        0123456789        0123456789
    In  : 1121112111        0010112111        0110112111        0110111111        0000001111        0110111101        0000000000
    Out : 1121112111        0010112111        0020112111        0010112111        0000001111        0010111111        0000000000
    source :                2                 3                 2                 2                 6                 2
    target :                1                 2                 6                 1                 8                 1
    branch :                [2]               []                [2]               [2,6]             [2,6]             [2,6]
    once selected :         {2:[1]}           {2:[1]}           {2:[1,6]}         {2:[1,6],6:[5]}   {2:[1,6],6:[5,8]}  "
    fail:                   [10321]           [10321]           [20321]           [10321,103265421]  "                 "
    path:                   10321->?          1032              10326             103265421->?      103268            1032687965421
                            backtrack                                             backtrack
  * requirements:
    ; many variables for keep track of states
    ; algorithm for backtracking

═════════════════════════════════════════════════

2. Fleury's Algorithm
* cut edge (bridge, isthmus): an edge whose removal disconnects a component of the graph
* condition for cut edge:
  ; If there is no backward edge from a subgraph to its ancestor (parent)
* steps:
  ; start at any vertex (node), remove edge to the next vertex
  ; choose any edge from this vertex, but not a cut edge
    then remove edge to the next vertex
  ; repeat the second step until the edge to the first start vertex
            (1)               (2)               (3)'              (3)               (4)'              (4)
            4 < 5             4 < 5             4 < 5             4 < 5             4 < 5             4 * 5
            v   ^             v   ^             v   ^             v   ^             v   *             *   *
        1 < 2 > 6 > 8     1 < 2 > 6 > 8     1 * 2 > 6 > 8     1 < 2 * 6 > 8     1 < 2 * 6 > 8     1 * 2 * 6 * 8
        *   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
        0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
        10                1032              10321->cut        10326             103265->cut       1032687965421
* requirements:
  ; an efficient algorithm to check connectivity for each edge deletion
  ; an efficient algorithm to find or avoid bridge (cut edge)

═════════════════════════════════════════════════

3. Eunice(?)'s algorithm (avoid bridge, not find)
  * StackExchagne: Determine if edge is a bridge in a graph
    https://math.stackexchange.com/questions/3965493/determine-if-edge-is-a-bridge-in-a-graph

  * Target node of bridge has min distance to sink node.
    So avoid an edge with min distance to sink node.

  * how to get distance dictionary?
    {'0':['3'],'1':['0'],'2':['1','6'],'3':['2'],'4':['2'],'5':['4'],'6':['5','8'],'7':['9'],'8':['7'],'9':['6']};
        target  source   distance  visited          dic
      --------------------------------------------------------------------------
      - 1       2        0         1                1:0
      - 2     ┌ 4        1         1,2              1:0,2:1
              └ 3        1         1,2              1:0,2:1
      ┌ 4       5        2         1,2,4            1:0,2:1,4:2
      └ 3       0        2         1,2,3            1:0,2:1,3:2
      ┌ 5       6        3         1,2,4,5          1:0.2:1,4:2,5:3
      └ 0      [1]stop   3         1,2,3,0          1:0,2:1,3:2,0:3                  (*) ┐
      - 6     ┌[2]stop                                                                   │
              └ 9        4         1,2,4,5,6        1:0.2:1,4:2,5:3,6:4                  │ merge
      - 9       7        5         1,2,4,5,6,9      1:0.2:1,4:2,5:3,6:4,9:5              │
      - 7       8        6         1,2,4,5,6,9,7    1:0.2:1,4:2,5:3,6:4,9:5,7:6          │
      - 8      [6]stop   7         1,2,4,5,6,9,7,8  1:0.2:1,4:2,5:3,6:4,9:5,7:6,8:7  (*) ┘

═════════════════════════════════════════════════

4. rotate circle at a node with multiple out-degrees
  * start from small circle, rotate it, then move to the next small circle
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

═════════════════════════════════════════════════

5. recursive algorithm (mysterious)
  * mysterious algorithm, clojure (by hadaarjan)
  * Sample algorithm by hadaarjan (Brilliant!)
    He (or She) focused on the recursive aspect of this probem.
    Structure:
    - Outside the recursion
      - a list for storing path
      - a starting node
    - base case (condition for ending recursion)
      - if selected key has no value
    - recursive case
    - visualization
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

═════════════════════════════════════════════════

6. DFS - Tree traversal (OOP)
7. Divide into small circles and put them back together.

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
    - How to notice bridges easily in Fleury's algorithm
      https://math.stackexchange.com/questions/2153841/how-to-notice-bridges-easily-in-fleurys-algorithm
  * Veblen's theorem
      https://en.wikipedia.org/wiki/Veblen%27s_theorem
  * bridge finding algorithms
    - Tarjans Strongly Connected Components algorithm | Graph Theory
      https://www.youtube.com/watch?v=TyWtx7q2D7Y
    - Kosaraju's algorithm
      https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm
      Strongly Connected Components Kosaraju's Algorithm Graph Algorithm
      Kosaraju algorithm introduction
      https://www.youtube.com/watch?v=Qdh6-a_2MxE     <---
      https://www.youtube.com/watch?v=Jb1XlDsr46o
      https://www.youtube.com/watch?v=HOOmetF56BI
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
  * Veblen's theorem
    https://en.wikipedia.org/wiki/Veblen%27s_theorem
    TODO: Is this related to multiple small cycles in a Eulerian cycle?
  * All Eulerian cycles
    - BEST theorem
  * javascript sytax
    - Object keys are always strings.
      https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Property_Accessors
      Using an integer as a key in an associative array in JavaScript
      https://stackoverflow.com/questions/2002923/using-an-integer-as-a-key-in-an-associative-array-in-javascript
*/

/************************************************
Hierholzer's algorithm
************************************************/

// get degrees
function get_degrees(graph) {
    /**
     * {str:[a]} -> {str:[int,int]}       // no tuple type in Javascript, object property is always string
     * returns In/Out-degrees
     * >>> get_degrees({ 0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6] })
     *     {'0':[1,1],'1':[1,1],'2':[2,2],'3':[1,1],'4':[1,1],'5':[1,1],'6':[2,2],'7':[1,1],'8':[1,1],'9':[1,1]}
     */
    let dic = {};
    for (const key in graph) {
        // count In-degrees
        for (const val of graph[key]) {
            if (!(val in dic)) dic[val] = [1, 0];
            else dic[val][0]++;
        }
        // count Out-degrees
        if (!(key in dic)) dic[key] = [0, graph[key].length];
        else dic[key][1] += graph[key].length;
    }
    return dic;
}

// check if graph is Eulerian
function check_eulerian_cycle(graph) {
    /**
     * {str:[a]} -> bool
     * returns true if graph is Eulerian cycle, else false (all nodes with even degrees)
     * >>> check_eulerian_cycle({0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]})
     *     true
     */
    const degrees = get_degrees(graph);
    for (const key in degrees) {
        const sum = degrees[key].reduce((a, b) => a + b, 0);
        if (sum % 2 != 0) return false;
    }
    return true;
}

// simple walk and construct graph
function walk(graph, node) {
    /**
     * ({str:[str]},str) -> [str]
     * returns a path or cycle by simple random walk
     */
    let cycle = [node];
    while (graph[node].length != 0) {
        const values = graph[node];
        const target = values[Math.floor(Math.random() * values.length)];
        cycle.push(target);
        graph[node] = values.filter(e => e != target);
        node = target;
    }
    return cycle;
}

// BA03F: Hierholzer's algorithm
function eulerian_hierholzer(graph) {
    /**
     * {a:[a]} -> [a]
     * returns a Eluerian cycle, destructive
     * >>> var graph = {'0':['3'],'1':['0'],'2':['1','6'],'3':['2'],'4':['2'],'5':['4'],'6':['5','8'],'7':['9'],'8':['7'],'9':['6']};
     * >>> eulerian_hierholzer(graph)
     *     ['6','5','4','2','1','0','3','2','6','8','7','9','6']
     */
    // first, get a cycle by random walk
    let nodes = Object.keys(graph);
    let source = nodes[Math.floor(Math.random() * nodes.length)];  // initial source node, randomly selected
    let cycle = walk(graph, source);
    // if there's remaining unvisited nodes, get another cycle, merge it to the previous cycle, again and again
    while (Object.keys(graph).map(e => graph[e].length === 0).every(Boolean) != true) {
        source = cycle.filter(e => graph[e].length != 0)[0];
        const new_cycle = walk(graph, source);
        const idx = cycle.indexOf(source);
        cycle = [...cycle.slice(0, idx + 1), ...new_cycle.slice(1), ...cycle.slice(idx + 1)];
    }
    return cycle;
}

/************************************************
Fleury’s Algorithm
************************************************/

function check_bridge(graph, source, target) {
    /**
     * ({str:[str]},str,str) -> bool
     * return true if a edge is bridge (cut edge, isthmus, ...)
     * TODO:
     */

}

function eulerian_fleury(graph) {
    /**
     * {a:[a]} -> [a]
     * returns a Eulerian cycle
     * TODO:
     */
}

/************************************************
Eunice's algorithm
https://math.stackexchange.com/questions/3965493/determine-if-edge-is-a-bridge-in-a-graph
************************************************/

// construct distance graph
function distance_graph(graph, dic, start, target, dist, visited) {
    /**
     * {str:[a]} -> {str:int}
     * returns a distances to sink node, backtracking and get distances of all nodes
     * from start node, construct distance dictionary
     * >>> var graph = {'0':['3'],'1':['0'],'2':['1','6'],'3':['2'],'4':['2'],'5':['4'],'6':['5','8'],'7':['9'],'8':['7'],'9':['6']};
     * >>> distance_graph(graph,{},'1','1',0,[]);
     *     {'0':3,'1':0,'2':1,'3':2,'4':2,'5':3,'6':4,'7':6,'8':7,'9':5}
     */
    // A. base case:
    if (visited.includes(target)) return dic;

    // B. recursive case:
    visited.push(target);
    dic[start] = 0;
    let sources = [];
    for (const key in graph) {
        if (graph[key].includes(target)) sources.push(key);
    }
    ++dist;
    for (const source of sources) {
        if (!(source in dic)) dic[source] = dist;
        else if (dic[source] > dist) dic[source] = dist;
    }
    for (const source of sources) {
        distance_graph(graph, dic, start, source, dist, visited);
    }
    return dic;
}

// BA03F: Eunice's algorithm
function eulerian_eunice(graph, start) {
    /**
     * {a:[a]} -> [a]
     * returns a Eulerian cycle, destructive
     * >>> var graph = {'0':['3'],'1':['0'],'2':['1','6'],'3':['2'],'4':['2'],'5':['4'],'6':['5','8'],'7':['9'],'8':['7'],'9':['6']};
     * >>> eulerian_eunice(graph, '1')
     *     ['1','0','3','2','6','8','7','9','6','5','4','2','1']
     */
    const distances = distance_graph(graph, {}, start, start, 0, []);
    let cycle = [start];
    let cur = start;
    while (Object.keys(graph).map(e => graph[e].length == 0).every(Boolean) != true) {
        if (graph[cur].length == 0) break;
        // A. if a node with one out-degree
        while (graph[cur].length == 1) {
            cur = graph[cur].pop();
            cycle.push(cur);
        }
        // B. if a node with no out-degree
        if (graph[cur].length == 0) break;
        // C. if a node with multiple out-degrees
        const candidates = graph[cur];
        let max_distance = -Infinity;
        let target;
        // avoid a node with min-distance to sink
        for (const candidate of candidates) {
            const dist = distances[candidate];
            if (dist > max_distance) {
                max_distance = dist;
                target = candidate;
            }
        }
        cycle.push(target);
        graph[cur] = graph[cur].filter(e => e !== target);
        cur = target;
    }
    return cycle;
}


/************************************************
rotate path at a node with multiple out-degrees
************************************************/

// BA03F: rotation
function eulerian_rotation(graph) {
    /**
     * {str:[str]} -> [str]
     * returns a Eulerian cycle by rotating cycle
     */
    const nodes = Object.keys(graph);
    start = nodes[Math.floor(Math.random() * nodes.length)];
    let cycle = [start];
    while (Object.keys(graph).map(e => graph[e].length == 0).every(Boolean) != true) {
        // if the first and last elements are the same
        if (cycle[0] === cycle[cycle.length - 1]) {
            // rotate the cycle for the last (==first) elem to have one or more values (key concept)
            while (!(cycle[0] in graph)) {
                cycle.push(cycle[1]);
                cycle = cycle.slice(1);
            }
        }
        // continue walking
        const source = cycle[cycle.length - 1];
        cycle.push(graph[source].pop());
        // if key has no value, delete key from graph
        if (graph[source].length === 0) delete graph[source];
    }
    return cycle;
}

/************************************************
recursive algorithm
************************************************/

// BA03F: closure
function eulerian_recursive(graph) {
    /**
     * {str:[str]} -> [str]
     * returns a Eulerian cycle by recursive algorithm
     */
    const nodes = Object.keys(graph);
    start = nodes[Math.floor(Math.random() * nodes.length)];
    let tour = [];
    function cycle(n) {
        while (graph[n].length > 0) {
            const next = graph[n].pop();
            cycle(next);
        }
        tour.push(n);
    }
    cycle(start);
    tour = tour.reverse();
    return tour;
}

/************************************************
random walk
************************************************/

// BA03F: random, by luck
function eulerian_random(graph) {
    /**
     * {str:[a]} -> [a]
     * returns a Eulerian cycle by luck (random walk), destructive
     */
    // backup original graph by deep-copy (Javascript deep-copy)
    const backup = JSON.parse(JSON.stringify(graph));
    // select a source node randomly, push it to cycle
    let nodes = Object.keys(graph);
    let source = nodes[Math.floor(Math.random() * nodes.length)];  // initial source node, randomly selected
    let cycle = [source];
    // loop until any key has non-empty value
    let trial = 0;
    while (Object.keys(graph).map(e => graph[e].length == 0).every(Boolean) != true) {
        const values = graph[source];
        let target = values[Math.floor(Math.random() * values.length)];  // target node, randomly selected
        cycle.push(target);                                              // add it to cycle
        graph[source] = graph[source].filter(item => item !== target);   // remove a node from graph
        source = target;                                                 // update focus
        // if cycle is closed prematurely
        ++trial;
        if (graph[source].length == 0 && Object.keys(graph).map(e => graph[e].length == 0).every(Boolean) != true) {
            graph = JSON.parse(JSON.stringify(backup));
            nodes = Object.keys(graph);
            source = nodes[Math.floor(Math.random() * nodes.length)];
            cycle = [source];
        }
    }
    console.log(`loop ${trial} times`);
    return cycle;
}

/************************************************
Break graph into minimum cycles, rejoin them into one large cycle
************************************************/

/************************************************
tree traversal, OOP
************************************************/

/************************************************
all Eulerian cycles
************************************************/

function is_disconnected_graphs(graph) {
    /**
     * {str:[str]} -> bool
     * returns true if a graph is disconnected
     * TODO:
     */
}

function is_simple_graph(graph) {
    /**
     * {str:[str]} -> bool
     * returns true if a graph is simple graph
     * TODO:
     */
}

function uvw_bypass(graph) {
    /**
     * {str:[str]} -> {str:[str]}
     * returns uvw-bypass graph
     * TODO:
     */
}

function all_eulerian_cycles(graph) {
    /**
     * {str:[str]} -> [[str]]
     * returns all Eulerian cycles
     * TODO:
     */
    let all_graph = new Set(graph);
}

/*
╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║  <-- simple graph
║         v <- a node with indegree larger than 1 in G        ║  <-- In-degree(node) > 1
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║  <-- uvw-bypass
║                 if NewGraph is connected                    ║  <-- connected/disconnected
║                     add NewGraph to AllGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║  <-- uvw-bypass to cycle
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝
*/

// helper fx
function read_graph(lines) {
    /**
     * [str] -> {int:[int]}  <-- This is wrong in Javascript! Object key must be a string type.
     * [str] -> {str:[str]}  <-- This is correct. Very strange.
     * returns a Eulerian graph
     */
    let graph = {};
    for (const line of lines) {
        const [source, val_str] = line.split(' -> ');
        let vals = [];
        for (const val of val_str.split(',')) {
            vals.push(val);
        }
        graph[parseInt(source)] = vals;
    }
    return graph;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03f.txt').toString().split("\n");
    const graph = read_graph(lines);

    // Hierholzer's Algorithm
    let startTime = performance.now()
    let cycle = eulerian_hierholzer(graph, '5');
    console.log(cycle.join('->'));
    console.log(`${performance.now() - startTime} milliseconds`)

    /*
    // Eunice's Algorithm
    startTime = performance.now()
    cycle = eulerian_eunice(graph, '5');
    console.log(cycle.join('->'));
    console.log(`${performance.now() - startTime} milliseconds`)
    */

    /*
    // ratation
    startTime = performance.now()
    cycle = eulerian_rotation(graph);
    console.log(cycle.join('->'));
    console.log(`${performance.now() - startTime} milliseconds`)
    */

    /*
    // recursive algorithm
    startTime = performance.now()
    cycle = eulerian_recursive(graph);
    console.log(cycle.join('->'));
    console.log(`${performance.now() - startTime} milliseconds`)
    */
}

// execute main function
main()