"""
Rosalind: BA3F (difficulty: 4/5)
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

[ Notice ]
- difficulty levels: 4/5
  1/5: solved with simple brute-force approach
  2/5: brute-force approach takes some time than better algorithms
  3/5: brute-force exceeds time limit (5 minutes)
  4/5: can't be solved with simple nested loops or brute-force algorithm
  5/5: requires deeper scientific knowledge

-------------------------------------------------

Plan 1.
            4 ← 5
            ↓   ↑                   [2], [6] have two directions
    0 → 3 →[2]→[6]→ 8 → 7 → 9       points for rolling back
     ↖   ↙     ↖       ↙
        1           ↖←↙

- required algorithms:
  (1) algorithm for detection: at node [2] and [6]
  (2) algorithm for correction: rolling back or something
      => ex) rotation (by Elmar Hinz)
      => ex) recursion (by hadaarjan)
  (3) algorithm for mamaging datastructure (list or dictionary)
- problem solving steps:
  a. problems
     * Is the given graph circular or not?
     * Complexity problem
       - What if a node has more than two edges? three, four, five ...
       - The fraph in this exercise is relatively simple.
         But...what if the given graph is many circles with a shared point?
         Can the function of solving this problem be applied to these complex graphs?
     * From which node should I start?
       - If graph is circular, it doens't matter at all.
       - If graph is not circular, it does really matter!
         How can I choose the correct starting node?
     * Python things
       - The list passed to argument of function will be modified.
         That means the original data can be changed after function call.
     * roll bak, how and when?
       - Where did the wrong selection take place?
       - How many selections were wrong?
       - How can I roll back to the time just before the wrong decision is made?
     * recursion, iteration, which is the better?
- simulation (writing down... I know... but implementation is not easy...)
  graph                                                              cycle
  ----------------------------------------------------------------------------------------
  {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [0]
  {      1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [03]
  {      1:[0],2:[1,6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [032]
  {      1:[0],2:[  6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [0321]      wrong selection
  {            2:[  6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [03210]     closed cycle!
  ----------------------------------------------------------------------------------------
  {      1:[0],2:[1,6],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [032]       roll back TODO: How?
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [0326]
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[  8],7:[9],8:[7],9:[6]}  [03265]     wrong selection
  {      1:[0],2:[1  ],      4:[2],      6:[  8],7:[9],8:[7],9:[6]}  [032654]
  {      1:[0],2:[1  ],                  6:[  8],7:[9],8:[7],9:[6]}  [0326542]
  {      1:[0],                          6:[  8],7:[9],8:[7],9:[6]}  [03265421]
  {                                      6:[  8],7:[9],8:[7],9:[6]}  [032654210] closed cycle!
  ----------------------------------------------------------------------------------------
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}  [0326]      roll back
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5  ],7:[9],8:[7],9:[6]}  [03268]
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5  ],7:[9],      9:[6]}  [032687]
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5  ],            9:[6]}  [0326879]
  {      1:[0],2:[1  ],      4:[2],5:[4],6:[5  ],                 }  [03268796]
  {      1:[0],2:[1  ],      4:[2],5:[4],                         }  [032687965]
  {      1:[0],2:[1  ],      4:[2],                               }  [0326879654]
  {      1:[0],2:[1  ],                                           }  [03268796542]
  {      1:[0],                                                   }  [032687965421]
  {                                                               }  [0326879654210] success

Plan 2.
- using loop and random choice
  (0,[3])          (0,3)
  (1,[0])          (1,0)
  (2,[1,6])        randomly select (2,1) or (2,6)
  (3,[2])          (3,2)
  (4,[2])     ->   (4,2)
  (5,[4])          (5,4)
  (6,[5,8])        randomly select (6,5) or (6,8)
  (7,[9])          (7,9)
  (8,[7])          (8,7)
  (9,[6])          (9,6)
- I implemented some functions, but they are so slow!
  How can I improve those functions?

Plan 3.
- using "CompositinoGraph_k(text)"
  In Composition Graph, the edges are isolated, meaning that no two edges share a node.
- Can I roll back with "Composition graph" rather than "De Bruijn graph"?
  (0,[3])          (0,3)
  (1,[0])          (1,0)
  (2,[1,6])        (2,1)
                   (2,6)
  (3,[2])          (3,2)
  (4,[2])     ->   (4,2)
  (5,[4])          (5,4)
  (6,[5,8])        (6,5)
                   (6,8)
  (7,[9])          (7,9)
  (8,[7])          (8,7)
  (9,[6])          (9,6)

  0->3->[2]->1->0                (x)  How can I roll back to node [2]
        [2]->[6]->5->4->2->1->0  (x)
             [6]->8->7->9->6->5->4->2->1->0

- simulation
  A. composition graph, all "unique" edges
  (0,3),(1,0),(2,1),(2,6),(3,2),(4,2),(5,4),(6,5),(6,8),(7,9),(8,7),(9,6)
  graph                                 cycle        next
  -------------------------------------------------------
  03 10 21 26 32 42 54 65 68 79 87 96   [0]             0
  ** 10 21 26 32 42 54 65 68 79 87 96   [03]            3
  -- 10 21 26 ** 42 54 65 68 79 87 96   [032]           2  No. candidates >= 2
  -- 10 ** 26 -- 42 54 65 68 79 87 96   [0321]          1
  -- ** -- 26 -- 42 54 65 68 79 87 96   [03210]         0  premature closing!
  -- 10 21 26 -- 42 54 65 68 79 87 96   [032]           2  roll back
  -- 10 21 ** -- 42 54 65 68 79 87 96   [0326]          6  No. candidates >= 2
  -- 10 21 -- -- 42 54 ** 68 79 87 96   [03265]         5
  -- 10 21 -- -- 42 ** -- 68 79 87 96   [032654]        4
  -- 10 21 -- -- ** -- -- 68 79 87 96   [0326542]       2
  -- 10 ** -- -- -- -- -- 68 79 87 96   [03265421]      1
  -- ** -- -- -- -- -- -- 68 79 87 96   [032654210]     1  premature closing!
  -- 10 21 -- -- 42 54 65 68 79 87 96   [0326]          6  roll back
  -- 10 21 -- -- 42 54 65 ** 79 87 96   [03268]         8
  -- 10 21 -- -- 42 54 65 -- 79 ** 96   [032687]        7
  -- 10 21 -- -- 42 54 65 -- ** -- 96   [0326879]       9
  -- 10 21 -- -- 42 54 65 -- -- -- **   [03268796]      6
  -- 10 21 -- -- 42 54 ** -- -- -- --   [032687965]     5
  -- 10 21 -- -- 42 ** -- -- -- -- --   [0326879654]    4
  -- 10 21 -- -- ** -- -- -- -- -- --   [03268796542]   2
  -- 10 ** -- -- -- -- -- -- -- -- --   [032687965421]  1
  -- ** -- -- -- -- -- -- -- -- -- --   [0326879654210] 0

  B. non-binary tree
      03         My ideas:
      |            a. get all paths, then select best (longest) path
      32           b. "rolling back"
     /  \           - remember the node with branching, store the state, ...
   21    26         - mysterious algorithms
   |    /  \           - rotation
  10  65    68         - recursion
      |     |       - construct non-binary tree -> Depth First Search (DFS)
      54    87      - dividing list into two: one value group, two or more values group
      |     |    Problems:
      42    79     - Getting all possible paths is inefficient.
      |     |      - List can be mutated every time it goes through for-loop.
      21    96     - If the list has changed, it is difficult to revert to the previous state.
      |     |      - Deep copying list is time-consuming process.
      10    65     - Creating a new complex data structure is not easy. (non-binary tree)
            |      - It requires at least two variable to remember state at a node with braching.
            54       (one for the index of the node, the other for index of the branch)
            |      - Composition graph can't be a dictionary, so it may be so slow during looping.
            42
            |
            21
            |
            10

C. splitting list, getting fragments (not good idea)
  [[(0,3),(1,0),(3,2),(4,2),(5,4),(7,9),(8,7),(9,6)], [(2,1),(2,6)], [(6,5),(6,8)]]
  [03 10 32 42 54 79 87 96]    vs.    [21 26]  [65 68]  <-- simplified version
  [03-32 10 54-42 79-96 87] - [21 65] => [03-32-21-10 54-42 79-96-65 87]
                              [21 65] => ...
  [03-32 10 54-42 79-96 87] - [21 68] => [03-32-21-10 54-42 79-96-68-87]
                              [21 68] => ...
  [03-32 10 54-42 79-96 87] - [26 65] => [03-32-26-65-54-42-10 79-96 87]
                              [26 65] => ...
  [03-32 10 54-42 79-96 87] - [26 68] => [03-32-26-68-87-79-96 54-42-10]
                              [26 68] => ...

Plan 4.
- adjacent matrix? any hidden idea?
      | 0  1  2  3  4  5  6  7  8  9
    -------------------------------
    0 |          *                       0 -> 3
    1 | *                                1 -> 0
    2 |    *              *              2 -> 1,6
    3 |       *                          3 -> 2
    4 |       *                          4 -> 2
    5 |             *                    5 -> 4
    6 |                *        *        6 -> 5,8
    7 |                            *     7 -> 9
    8 |                      *           8 -> 7
    9 |                   *              9 -> 6

- simpler case
  {0:[1],1:[2],2:[3,0],3:[4],4:[2]}      0 -> 1
                 ↙ ↖                   1 -> 2
      0 → 1 →[2]→ 3 → 4                  2 -> 0,3
        ↖ ↙                            3 -> 4
                                         4 -> 2
      | 0 1 2 3 4     idx   01234
    ------------- ->  -------------
    0 |   *           0   [[01000]
    1 |     *         1    [00100]
    2 | *     *       2    [10010]
    3 |         *     3    [00001]
    4 |     *         4    [00100]]

    cycle   next (just referencing list, not mutating)
    ----------------------------------------------------------
    0       1
    01      2
    012     0 already exists, and has just one target -> next 3
    0123    4
    01234   2 already exists, but has two targets -> next 2
    012342  0 already exists, next 3 already exists -> end!

References:
- Eulerian path
  https://en.wikipedia.org/wiki/Eulerian_path
- Fleury's Algorithm
  https://replit.com/@PascalTherese/Fleurys-Algorithm
- Hierholzer’s Algorithm for directed graph
  https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
- Delete an element from a dictionary
  https://stackoverflow.com/questions/5844672/delete-an-element-from-a-dictionary
- Python: Checking if a 'Dictionary' is empty doesn't seem to work
  https://stackoverflow.com/questions/23177439/python-checking-if-a-dictionary-is-empty-doesnt-seem-to-work
- What is the difference between the random.choices() and random.sample() functions?
  https://stackoverflow.com/questions/59763933/what-is-the-difference-between-the-random-choices-and-random-sample-function
- Deep copy of a dict in python
  https://stackoverflow.com/questions/5105517/deep-copy-of-a-dict-in-python
- How do I pass a variable by reference?
  https://stackoverflow.com/questions/986006/how-do-i-pass-a-variable-by-reference
- How do I write a function with output parameters (call by reference)?
  https://docs.python.org/3/faq/programming.html#how-do-i-write-a-function-with-output-parameters-call-by-reference
  Remember that arguments are passed by assignment in Python.
  Since assignment just creates references to objects,
  there’s no alias between an argument name in the caller and callee,
  and so no call-by-reference per se.
  You can achieve the desired effect in a number of ways.
- What's the idiomatic syntax for prepending to a short python list?
  https://stackoverflow.com/questions/8537916/whats-the-idiomatic-syntax-for-prepending-to-a-short-python-list
- Why aren't python nested functions called closures?
  https://stackoverflow.com/questions/4020419/why-arent-python-nested-functions-called-closures
- Python Inner Functions: What Are They Good For?
  https://realpython.com/inner-functions-what-are-they-good-for/
- How can I implement a tree in Python?
  https://stackoverflow.com/questions/2358045/how-can-i-implement-a-tree-in-python
- Non-Binary Tree Data Structure in Python
  https://stackoverflow.com/questions/60579330/non-binary-tree-data-structure-in-python
- how to find a path in non-binary tree python
  https://stackoverflow.com/questions/67164207/how-to-find-a-path-in-non-binary-tree-python
- PEP 498 – Literal String Interpolation
  https://peps.python.org/pep-0498/
"""
#!/usr/bin/env python
import time
import random


def read_de_bruijn_graph_as_list(lines):
    """
    De Bruijn Graph as list
    >>> read_de_bruijn_graph_as_list(lines)
        [(0,[3]),(1,[0]),(2,[1,6]),(3,[2]),(4,[2]),(5,[4]),(6,[5,8]),(7,[9]),(8,[7]),(9,[6])]
    """
    adjlist = []
    for line in lines:
        start, temp = [line.strip() for line in line.split(' -> ')]
        ends = [int(elem) for elem in temp.split(',')]
        adjlist.append((int(start), ends))
    return adjlist


def read_de_bruijn_graph_as_dict(lines):
    """
    De Bruijn Graph as dictionary
    >>> read_de_bruijn_graph_as_dict(lines)
        {0: [3],1: [0],2: [1,6],3: [2],4: [2],5: [4],6: [5,8],7: [9],8: [7],9: [6]}
    """
    adjlist = {}
    for line in lines:
        start, temp = [line.strip() for line in line.split(' -> ')]
        ends = [int(elem) for elem in temp.split(',')]
        adjlist[int(start)] = ends
    return adjlist


def read_debruijn_to_composition_as_list(lines):
    """
    CompositionGraph
    >>> read_debruijn_to_composition_as_list(lines)
        [(0,3),(1,0),(2,1),(2,6),(3,2),(4,2),(5,4),(6,5),(6,8),(7,9),(8,7),(9,6)]
    """
    adjlist = []
    for line in lines:
        start, temp = [line.strip() for line in line.split(' -> ')]
        ends = [int(elem) for elem in temp.split(',')]
        for end in ends:
            adjlist.append((int(start), int(end)))
    return adjlist


def eulerian_01_incomplete(adjdict):
    # choose a key in adjacent dictionary
    start = min(adjdict.keys())
    cycle = [start]
    # loop until cycle is correct
    while adjdict:
        # if the cycle ends prematurely
        if cycle[0] == cycle[-1]:
            # do something for rolling back
            # TODO: Implement this mysterious part
            pass
        # select the next node
        lst = cycle[-1]
        cycle.append(adjdict[lst])
        # if key has no value
        if not adjdict[lst]:
            del adjdict[lst]
    return cycle


def eulerian_pop_incomplete(adjdict):
    """ wrong, TODO: Fix this
    This function sometimes returns correct cycle, but not always.
    """
    start = min(adjdict.keys())
    cycle = [start]
    while adjdict:
        if start not in adjdict:
            return cycle
        end = adjdict[start].pop()  # <--- always extract the last one
        cycle.append(end)
        if adjdict[start] == []:
            del adjdict[start]
        start = end
    return cycle


def eulerian_rand_ver01_failed(adjdict):
    """ wrong, TODO: Fix this
    This function sometimes returns correct cycle, but not always.
    """
    start = min(adjdict.keys())
    cycle = [start]
    while adjdict:
        if start not in adjdict:
            break
        if len(cycle) == 13:
            break
        end = random.sample(adjdict[start], k=1)[0]  # <-- randomly select path
        adjdict[start].remove(end)
        cycle.append(end)
        if adjdict[start] == []:
            del adjdict[start]
        start = end
    return cycle


def eulerian_rand_ver02_slow(adjdict):
    """ it works, but too slow TODO: Fix this """
    start = min(adjdict.keys())
    cycle = [start]
    while True:
        copied = {key: value[:] for key, value in adjdict.items()}
        temp_cycle = [start]
        while copied:
            if start not in copied:
                break
            end = random.sample(copied[start], k=1)[0]
            copied[start].remove(end)
            temp_cycle.append(end)
            if copied[start] == []:
                del copied[start]
            start = end
        cycle = temp_cycle
        if copied == {}:
            return cycle


def eulerian_rand_ver03_slow(adjdict):
    """ it works, but too slow TODO: Fix this """
    while True:
        copied = {key: value[:] for key, value in adjdict.items()}
        cycle = eulerian_rand_ver01(copied)
        if not copied:
            return cycle


def eulerian_recur01_incomplete(adjdict, start, output):
    """
    just for testing
    This function only work when
    every key in de Bruijn graph (as dictionary) has just one value.
    0->1->2->3->4->0
    >>> graph = {0:1, 1:2, 2:3, 3:4, 4:0}
    >>> cycle = eulerian_recur01(graph, 0, [0])
        [0, 1, 2, 3, 4, 0]
    """
    if not adjdict:
        return output
    lst = output[-1]
    output.append(adjdict[lst])
    del adjdict[lst]
    return eulerian_recur01(adjdict, start, output)


def eulerian_recur02_todo(adjdict, start, output):
    """
    just for testing
      ↙  ↖
    4 → 0 → 1 → 2 → 3
              ↖  ↙
    >>> adjdict = {0:[1], 1:[2,4], 2:[3], 3:[1], 4:[0]}
    >>> cycle = eulerian_recur02(adjdict, 0, [0])
    """
    pass


def main():
    f = open('/home/wsl/rosalind/data/ba03f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    gdict = read_de_bruijn_graph_as_dict(lines)
    samplegraph = {0:[1], 1:[2], 2:[3,0], 3:[4], 4:[2]}

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()