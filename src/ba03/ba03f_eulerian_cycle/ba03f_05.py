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

═════════════════════════════════════════════════

- Sample algorithm by Elmar Hinz (Very clever!)
  Example A.
                 ↙ ↖
      0 → 1 →[2]→ 3 → 4
        ↖ ↙
  graph = {0:[1], 1:[2], 2:[3,0], 3:[4], 4:[2]}
  how rotation works when cycle closes:
  cycle               ramining graph
  --------------------------------------------
  0->1->2->0          {2:[3],3:[4],4:[2]}
     1->2->0          {2:[3],3:[4],4:[2]}  delete the 1st node
     1->2->0->1       {2:[3],3:[4],4:[2]}  append the 1st node to the end (no further path)
        2->0->1       {2:[3],3:[4],4:[2]}  delete the 1st node
        2->0->1->2    {2:[3],3:[4],4:[2]}  append the 1st node to the end (has further path)

  Example B.
            4 ← 5
            ↓   ↑
    0 → 3 →[2]→[6]→ 8 → 7 → 9
     ↖   ↙     ↖       ↙
        1           ↖←↙

  cycle
  --------------------------------------------
  [0]               {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [03]              {      1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [032]             {      1:[0],2:[1,6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [0321]            {      1:[0],2:[  6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [03210]           {            2:[  6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [-3210]           the first rotation starts here
  [32103]           ... delete the first, the next first to the end ...
  [-2103]
  [21032]           the first rotation ends here
  [210326]          {                          4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [210326]          {                          4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
  [2103265]         {                          4:[2],5:[4],6:[  8],7:[9],8:[7],9:[6]}
  [21032654]        {                          4:[2],      6:[  8],7:[9],8:[7],9:[6]}
  [210326542]       {                                      6:[  8],7:[9],8:[7],9:[6]}
  [-10326542]       the second rotation starts here
  [103265421]       ... delete the first, the next first to the end ...
  [-03265421]
  [032654210]
  [326542103]
  [265421032]
  [654210326]
  [542103265]
  [421032654]
  [210326542]
  [103265421]
  [032654210]
  [326542103]
  [265421032]      the second rotation ends here
  [654210326]      {                                      6:[  8],7:[9],8:[7],9:[6]}
  [6542103268]     {                                              7:[9],8:[7],9:[6]}
  [65421032687]    {                                              7:[9],      9:[6]}
  [654210326879]   {                                                          9:[6]}
  [6542103268796]  {                                                               }
'''

def eulerianCycle_Hinz(graph):
    ''' by Elmar Hinz '''
    # init cycle with any node
    cycle = [min(graph.keys())]
    # run until all edges are moved from graph to cycle
    while len(graph) > 0:
        # whenever cycle closes, rotate to a node with remaining targets
        if cycle[0] == cycle[-1]:
            while not cycle[0] in graph:
                cycle.pop(0)
                cycle.append(cycle[0])
        # the last node of cycle is the new source
        source = cycle[-1]
        # move one target at a time from graph to the end of cycle
        cycle.append(graph[source].pop())
        # clean up empty dict entries of graph
        if len(graph[source]) == 0: del graph[source]
    return cycle