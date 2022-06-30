"""
Rosalind: BA6C
Compute the 2-Break Distance Between a Pair of Genomes

2-Break Distance Problem
Find the 2-break distance between two genomes.

Given: Two genomes with circular chromosomes on the same set of synteny blocks.

Return: The 2-break distance between these two genomes.

Sample Dataset
(+1 +2 +3 +4 +5 +6)
(+1 -3 -6 -5)(+2 -4)

Sample Output
3

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> Cycle to Chromosome (BA6G)
      -> Genome to ColoredEdges (BA6H)
      -> ColoredEdges to Genome (BA6I)
      -> 2-Break on ColoredEdges (BA6J)
      -> PREV: 2-Break on Genome (BA6K)
      ↓
    * HERE: 2-Break distance (BA6C)
      ↓
    * NEXT: 2-Break sorting (BA6D)

Info.
  * This quiz requires lots of background knowledge.
    I had to draw a lot of circles on the paper to solve this quiz.
    This quiz is not easy at all.

  * genome to graph - fig 6.9
    - Circular chromosome with n elements can be written in 2n different ways.
      circle 1            graph representations
      -----------------------------------------
      (+a -b -c +d)       (+a -b -c +d)
                          (-b -c +d +a)
                          (-c +d +a -b)
                          (+d +a -b -c)
                          (-a -d +c +b)
                          (-d +c +b -a)
                          (+c +b -a -d)
                          (+b -a -d +c)

  * BreakPointGraph(P,Q) - fig 6.14~16
    - How to draw:
      a. Let's pick one Chromosome to be the reference, and let's say it's chromosome P.
      b. Draw a circle.
      c. Mark the position and direction of the synteny block with an arc on the circle.
         The arc connecting the synteny blocks is painted Blue (colored edges).
      d. Repeat the process from (a) to (c) to draw another new chromosome Q.
         The arc connecting the synteny blocks is painted Red (colored edges).
      e. Then, overlap two circles.
         When two circles are overlapped,
         Blue and Red arcs are connected to each other
         because the directions and positions of the synteny blocks of the two circles are drawn the same.
         That is "BreakPointGraph(P,Q)",
         and circles that are made of red and blue lines are the "Cycles(P,Q)".
         Cycles of length 2 are "trivial cycles".
    - Colored edges of P and Q is BreakPointGraph(P,Q).
      That means that BreakPointGraph(P,Q) is a union of two sets.
      BreakPointGraph(P,Q) = ColoredEdges(P) + ColoredEdges(Q)
    - "Cycles(P,Q)" means circular graph (circle) in BreakPointGraph(P,Q).
    - "Block(P,Q)": the number of synteny block

  * 2-breaks - fig 6.10
    - definition: general operation on the genome graph
                  in which two red edges are replaced
                  with two new red edges on the same four nodes
    - the result of 2-breaks
                      P(+a -b -c +d)     Q(+a -b -d +c)
      (1) reversal     (+a -b -c +d)      (+a -b -d +c)
      (2) fission      (+a -b -c +d)      (+a -b)(-c +d)
      (3) fusion       (+a -b)(-c +d)     (+a -b -c +d)
      (4) translocation
                      P(-a +b +c -d)     Q(+e +f -g +h)
                 -->  T(+e +b +c -d -a +f -g +h)

           -a      +b      +c      -d
         P <──*─x─*──>*───*──>*───*<──   two chromosomes
         Q ──>*─x─*──>*───*<──*───*──>
           +e      +f      -g      +h

           into

           -a      +f      -g      +h
           <──*───*──>*───*<──*───*──>   one chromosomes
           ──>*───*──>*───*──>*───*<──
           +e      +b      +c      -d

  * 2-breaks distance
    - definition: the number of operations in a shortest sequence of 2-breaks
                  transforming P into Q as the 2-break distance
                  between P and Q, denoted d(P,Q)
    - d(P,Q) = Block(P,Q) - Cycles(P,Q)

  * theorems
    (1) Breakpoint theorem:
        d_rev(P) is greater than or equal to BreakPoint(P)/2
    (2) Cycle theorem:
        for genomes P and Q, any 2-break applied to P can increase Cycles(P,Q) by at most 1
    (3) 2-break distance theorem:
        the 2-break distance between genomes P and Q is equal to Block(P,Q) - Cycles(P,Q)
        d(P,Q) = Block(P,Q) - Cycles(P,Q)

  * steps:
    - get the number of synteny blocks from chromosome P, Q --> Blocks(P,Q)
    - get a set of colored edges from genome P
    - get a set of colored edges from genome Q
    - combine two sets into one --> BreakPointGraph(P,Q)
    - get the number of cycles in the combined set --> Cycles(P,Q)  <-- key component in this quiz
    - d(P,Q) = Blocks(P,Q) - Cycles(P,Q)

Plan 1. (Failed)
  * using 1:1 dictionary
  * sample dataset
    Blocks(P,Q) = 6
    P: (+1 +2 +3 +4 +5 +6)  -> {(10,11),(12,1),(2,3),(6,7),(4,5),(8,9)}
    Q: (+1 -3 -6 -5)(+2 -4) -> {(11,10),(7,3),(2,6),(4,8),(9,1),(5,12)}
    BrealPointGraph(P,Q): {(10,11),(12,1),(11,10),(7,3),(2,3),(6,7),(4,5),(8,9),(2,6),(4,8),(9,1),(5,12)}
    Cycles(P,Q): 10-11-10, 1-9-8-4-5-12-1, 6-7-3-2-6

  * Cycles(P,Q)
    (1) make the first element unique integer (key)
    (2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1)
     2:6   3:2   4:5   8:4   5:12   6:7   7:3   9:8   1:9   10:11   11:10   12:1  <-- unique keys
     1:9,2:6,3:2,4:5,5:12,6:7,7:3,8:4,9:8,10:11,11:10,12:1

  * problem:
    a. set --> dictionary, conversion problem (undirected graph!)
    {(2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1)}  <-- sorted
    {1:9,9:8,8:4,4:5,5:12,12:1,10:11,11:10,7:3,3:2,2:6,6:7}  --> It works.

    {(2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1)}
    {10:11,12:1,11:10,7:3,2:3,4:5,6:2,8:4,9:1,5:12}  --> It down't works.

Plan 2. (Stopped)
  * adding edges with reverse directions
  * undirected edges --> directed edges, by adding edges with reversed directions
    {(2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1),
     (6,2),(3,2),(5,4),(8,4),(12,5),(7,6),(3,7),(9,8),(1,9),(11,10),(10,11),(1,12)}

  * select one edge, do something, then remove three edges
    {(2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1),
     (6,2),(3,2),(5,4),(8,4),(12,5),(7,6),(3,7),(9,8),(1,9),(11,10),(10,11),(1,12)}
    selected: (2,6)
    find the next edge <-- requires another loop! problem
    removed: (6,2), (6,7), (7,6)

Plan 3. (Solved)
  * undirected --> bidirectional graph as dictionary
    {(2,6),(2,3),(4,5),(4,8),(5,12),(6,7),(7,3),(8,9),(9,1),(10,11),(11,10),(12,1)}
      2:6   2:3   4:5   4:8   5:12   6:7   7:3   8:9   9:1   10:11   11:10   12:1
      6:2   3:2   5:4   8:4   12:5   7:6   3:7   9:8   1:9   11:10   10:11   1:12
    {1:[9,12],2:[3,6],3:[2,7],4:[5,8],5:[4,12],6:[2,7],7:[3,6],8:[4,9],9:[1,8],10:[11,11],11:[10,10],12:[1,5]}
    --> solved, but not efficient

═════════════════════════════════════════════════

References:
- Difference between union() and update() in sets, and others?
  https://stackoverflow.com/questions/13905640/difference-between-union-and-update-in-sets-and-others
- Pythonic way to access arbitrary element from dictionary [duplicate]
  https://stackoverflow.com/questions/10593651/pythonic-way-to-access-arbitrary-element-from-dictionary
- Access an arbitrary element in a dictionary in Python - popitem()
  https://stackoverflow.com/questions/3097866/access-an-arbitrary-element-in-a-dictionary-in-python
- How can I verify if one list is a subset of another?
  https://stackoverflow.com/questions/16579085/how-can-i-verify-if-one-list-is-a-subset-of-another
- One liner to determine if dictionary values are all empty lists or not
  https://stackoverflow.com/questions/5889611/one-liner-to-determine-if-dictionary-values-are-all-empty-lists-or-not
- Python Closures
  https://www.programiz.com/python-programming/closure
"""
#!/usr/bin/env python
import time


def allpaths(graph, start, sink, path=[]):
    """
    ({a:[a]},a,[a]) -> [[a]]
    returns all path from 'start' node to 'sink' node (BA5D)
    TODO: Is this useful in this exercise?
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


def chromosome_to_cycle(chromosome):
    """
    [int] -> [int]
    algorithm in textbook
    return a list of integers
    """
    nodes = []
    for j in range(len(chromosome)):
        i = chromosome[j]
        if i > 0:
            nodes.append(2*i - 1)
            nodes.append(2*i)
        else:
            nodes.append(-2*i)
            nodes.append(-2*i - 1)
    return nodes


def cycle_to_chromosome(cycle):
    """
    [int] -> [int]
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])
        [1,-2,-3,4]
    >>> cycle_to_chromosome([2,4,3,6,5,1,2])
        [1,-2,-3]
    >>> cycle_to_chromosome([2,3,4,8,7,6,5,1,2])
        [2,-4,-3,1]

    """
    if cycle[0] == cycle[-1]:
        cycle = cycle[:-1]
    while not (cycle[0]%2 == 1 and cycle[1]%2 == 0 and (cycle[0]-1 == cycle[1] or cycle[0]+1 == cycle[1])):
        cycle = cycle[1:] + [cycle[0]]
    chromosome = []
    for i in range(0, len(cycle), 2):
        if cycle[i] < cycle[i + 1]:
            chromosome.append(cycle[i + 1] // 2)
        else:
            chromosome.append(-cycle[i] // 2)
    return chromosome


def genome_to_colored(genome):
    """
    [[int]] -> {(int,int)}
    returns a list of list of integers
    >>> genome_to_colored([[1,-2,-3],[4,5,-6]])
        {(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)}
    """
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


def two_breakpoint_graph(genome_p, genome_q):
    """
    ([[int]],[[int]]) -> {(int,int)}
    implementation of 'BreakPointGraph(P,Q)' function in textbook
    Is this just union of two sets? No.
    '2-BreakOnGenomeGraph' function in BA6J can create new red-colored edges already in blue-colored edges.
    blue = {(2,3),(7,6),(5,1),(4,8)}
    red = {(2,4),(3,6),(5,7),(8,1)}
    two = 2-BreakOnGenomeGraph(red,2,4,3,6) --> {(2,3),(4,6),(5,7),(8,1)}
    blue.union(red) --> {(2,3),(7,6),(5,1),(4,8),(3,6),(5,7),(8,1)}   <-- one (2,3) was removed
    It is because of trivial cycle, and it can be a problem.
    {(2,3)} is nothing at all, but {(2,3),(3,2)} is a trivial cycle.
    That a big difference.
    """
    colored_p = genome_to_colored(genome_p)
    colored_q = genome_to_colored(genome_q)
    colored = colored_p.union(colored_q)
    return colored


def num_cycles_wrong1(graph):
    """
    {(a,a)} -> int
    implementation of 'Cycles(P,Q)' function in textbook
    TODO: It works for sample dataset, but TOO SLOW. Sometimes incorrect. Make it better.
    """
    cycles = []
    cur = graph.pop()
    cycle = [cur[0], cur[1]]
    while graph:
        if cycle[0] == cycle[-1]:
            cycles.append(cycle)
            cur = graph.pop()
            cycle = [cur[0], cur[1]]
        for edge in graph:
            # append
            if cycle[-1] == edge[0]:
                cur = edge
                cycle = cycle + [cur[1]]
                break
            elif cycle[-1] == edge[1]:
                cur = edge
                cycle = cycle + [cur[0]]
                break
            elif cycle[0] == edge[0]:
                cur = edge
                cycle = [cur[1]] + cycle
                break
            elif cycle[0] == edge[1]:
                cur = edge
                cycle = [cur[0]] + cycle
                break
        graph.remove(cur)
        if not graph:
            cycles.append(cycle)
    return len(cycles)


def num_cycles_wrong2(edges):
    """
    {(a,a)} -> int
    returns the number of cycles in a set of edges
    """
    # edges to graph with unique keys
    graph = {}
    # edges = sorted(edges)
    for edge in edges:
        if edge[0] not in graph:
            graph[edge[0]] = edge[1]
        else:
            graph[edge[1]] = edge[0]
    print(graph)
    # get the number of cycles
    count = 1
    within_cycle = False
    next_node = 0
    while graph:
        if not within_cycle:
            _, next_node = graph.popitem()
            within_cycle = True
        if next_node not in graph:
            count += 1
            within_cycle= False
            continue
        next_node = graph.pop(next_node)
    return count


def num_cycles(edges):
    """
    {(a,a)} -> {a:[a]} -> int
    returns the number of cycles from a undirected graph (a set of tuples)
    >>> num_cycles_wrong3({(1,9),(2,6),(3,2),(4,5),(5,12),(6,7),(7,3),(8,4),(9,8),(10,11),(11,10),(12,1)})
        3
    TODO: Too slow! Make this fx better.
    """
    # a set of undirected edges -> a dictionary of bidirectional edges
    graph = {}
    for edge in edges:
        if edge[0] not in graph:
            graph[edge[0]] = [edge[1]]
            if edge[1] not in graph:
                graph[edge[1]] = [edge[0]]
            else:
                graph[edge[1]].append(edge[0])
        else:
            graph[edge[0]].append(edge[1])
            if edge[1] not in graph:
                graph[edge[1]] = [edge[0]]
            else:
                graph[edge[1]].append(edge[0])
    # the number of cycles
    count = 1
    in_cycle = False
    next_key = 0
    while any([graph[i] != [] for i in graph]):
        # start a new cycle
        if not in_cycle:
            cur_key = [k for k,v in graph.items() if v][0]
            cur_val = graph[cur_key].pop()
            next_key = cur_val
            if not graph[cur_key]:
                del graph[cur_key]
            if cur_val in graph:
                graph[cur_val].remove(cur_key)
                if not graph[cur_val]:
                    del graph[cur_val]
            else:
                del graph[cur_val]
            in_cycle = True
        # end of a cycle
        if next_key not in graph:
            in_cycle = False
            count += 1
            continue
        # continue cycle
        popped = graph[next_key].pop()
        if not graph[next_key]:
            del graph[next_key]
        if popped in graph:
            graph[popped].remove(next_key)
            if not graph[popped]:
                del graph[popped]
        else:
            del graph[popped]
        next_key = popped
    return count


def two_break_distance(genome_p, genome_q):
    '''
    ([[int]],[[int]]) -> int
    returns 2-break distance (d(P,Q))
    '''
    colored_p = genome_to_colored(genome_p)
    colored_q = genome_to_colored(genome_q)
    colored = colored_p.union(colored_q)
    return sum(map(len, genome_p)) - num_cycles(colored)


def read_genome(line):
    """
    str -> [[int]]
    returns a genome from a string line (helper function)
    >>> read_genome("(+1 -3 -6 -5)(+2 -4)")
        [[1,-3,-6,-5],[2,-4]]
    >>> read_genome("(+1 +2 +3 +4 +5 +6)")
        [[1,2,3,4,5,6]]
    """
    genome = []
    chromosome = []
    elem = ''
    for ch in line:
        if ch == '(':
            continue
        elif ch not in "(), ":
            elem += ch
        elif ch == ')':
            chromosome.append(int(elem))
            genome.append(chromosome)
            chromosome = []
            elem = ''
        else:
            chromosome.append(int(elem))
            elem = ''
    return genome


def main():
    f = open('/home/wsl/rosalind/data/ba06c.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    genome_p = read_genome(lines[0])
    genome_q = read_genome(lines[1])
    f.close()

    sample_genome_p = [[+1,+2,+3,+4,+5,+6]]
    sample_genome_q = [[+1,-3,-6,-5],[+2,-4]]

    start_time = time.time()
    print(two_break_distance(genome_p,genome_q))
    print("--- %s seconds ---" % (time.time() - start_time))  # ~ 40 seconds!! too slow


if __name__ == "__main__":
    main()
