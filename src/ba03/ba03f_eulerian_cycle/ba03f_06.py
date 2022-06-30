'''
Rosalind: BA3F
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

    EULERIANCYCLE(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle' by traversing Cycle (starting at newStart) and then randomly walking
            Cycle <- Cycle'
        return Cycle

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

-------------------------------------------------

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

    def eulerianCycle_hadaarjan(gdict):
        start = min(gdict.keys())
        tour = []
        def cycle(n):
            print("A: ", n, tour, gdict)
            while len(gdict[n]) > 0:
                next = gdict[n].pop()
                print("B: ", n, next, gdict)
                cycle(next)
            print("C: ", n, tour, gdict)
            tour.append(n)
        cycle(start)
        tour = tour[::-1]
        print(tour)
        return tour

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
'''

def eulerianCycle_hadaarjan(gdict):
    ''' by hadaarjan '''
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