"""
Rosalind: BA6D
Find a Shortest Transformation of One Genome into Another by 2-Breaks

2-Break Sorting Problem
Find a shortest transformation of one genome into another by 2-breaks.

Given:
Two genomes with circular chromosomes on the same set of synteny blocks.

Return:
The sequence of genomes resulting from applying a shortest sequence of 2-breaks
transforming one genome into the other.

Sample Dataset
(+1 -2 -3 +4)
(+1 +2 -4 -3)

Sample Output
(+1 -2 -3 +4)
(+1 -2 -3)(+4)
(+1 -2 -4 -3)
(+1 +2 -4 -3)

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
      -> 2-Break on Genome (BA6K)
      ↓
    * PREV: 2-Break distance (BA6C)
      ↓
    * HERE: 2-Break sorting (BA6D)
      ↓
    * NEXT: Shared k-mers (BA6E)

Info.
  * This quiz lacks information. Give me some more info!
    This quiz can't be solved without textbook.
    Even with textbook, this quiz is not easy at all.
    - What is the "2-break sorting"?
      I understand the concept of 2-break, but what is the meaning of "2-break sorting"?
    - the format of sample output?
    - the order of genomes in chromosome?
    - the order of chromosome in genome?
    - the first genome of a chromosome?
    - the first chromosome of a genome?
    - sign of the first genome of a chromosome?
    - break two red-colored edges in nontrivial cycle, increase the number of trivial cycles...
      Is that all? Do I need to learn more knowledge academically?

  * algorithm in textbook
  ╔════════════════════════════════════════════════════════════════════════════════════╗
  ║ SHORTESTREARRANGEMENTSCENARIO(P, Q)                                                ║
  ║     output P                                                                       ║
  ║     RedEdges <- COLOREDEDGES(P)                                                    ║
  ║     BlueEdges <- COLOREDEDGES(Q)                                                   ║
  ║     BreakpointGraph <- the graph formed by RedEdges and BlueEdges                  ║
  ║     while BreakpointGraph has a non-trivial cycle Cycle                            ║
  ║         (j, i') <- an arbitrary edge from BlueEdges in a nontrivial red-blue cycle ║
  ║         (i, j) <- an edge from RedEdges originating at node j                      ║
  ║         (i', j') <- an edge from RedEdges originating at node i'                   ║
  ║         RedEdges <- RedEdges with edges (i, j) and (i', j') removed                ║
  ║         RedEdges <- RedEdges with edges (j, i') and (j', i) added                  ║
  ║         BreakpointGraph <- the graph formed by RedEdges and BlueEdges              ║
  ║         P <- 2-BREAKONGENOME(P, i, i', j, j')   <-- TODO: Is this correct?         ║
  ║     output P                                                                       ║
  ╚════════════════════════════════════════════════════════════════════════════════════╝

  * review oof theorems in chapter 06
    Theorems:
    (1) Breakpoint theorem:
        d_rev(P) is greater than or equal to BreakPoint(P)/2
    (2) Cycle theorem:
        for genomes P and Q, any 2-break applied to P can increase Cycles(P,Q) by at most 1
    (3) 2-break distance theorem:
        the 2-break distance between genomes P and Q is equal to Block(P,Q) - Cycles(P,Q)
        d(P,Q) = Block(P,Q) - Cycles(P,Q)

  * making all cycles trivial == make P equal to Q
    break (remove) two red-colored edges, then add new two red-colored edges

  * How can I successfuly break two colored edges in P
    and increase the number of trivial cycles by 1?
    ; Select two red edges with one blue edge between them.
                       red     blue     red
      before:      ip────────i────────j────────jp

                       blue                 red
      after :       i════════j       ---ip────────jp---
                       red
                   In this way, I can increase the number of trivial cycles.

  * cycle vs. geneome
    i.  cycles --> I need 'BreakPointGraph(P,Q)' which is a set of combined blue-colored and red-colored edges.
    ii. genome --> I need only a set of red-colored edges.
    So I have to do two separated processes at the same time in the 'SHORTESTREARRANGEMENTSCENARIO'.
       (1) increasing the number of trivial cycles
           == merging two sets (red/blue-colored edges) into 'BreakPointGraph(P,Q)'
              then, get nontrivial cycles from the 'BreakPointGraph(P,Q)'
              Red-colored edges should not be modified. (If so, it can make errors in process (2).)
       (2) changing P to P' again and again until there's no nontrivial cycle
           == applying 2-BreakOnGenomeGraph to a set of red-colored edges

  * sample dateset
    P: (+1 -2 -3 +4)  --> {(2,4),(8,1),(5,7),(3,6)}
    Q: (+1 +2 -4 -3)  --> {(2,3),(7,6),(5,1),(4,8)}

    genome  P           red-colored edges              cycles                diagram                             num of nontrivial/trivial cycles
    ------------------------------------------------------------------------------------------------------------------------------
    (+1 -2 -3 +4)   --> P0: {(2,4),(3,6),(8,1),(5,7)}  1-2-4-3-6-5-7-8-1     1──>2───4<──3───6<──5───7──>8───1   1 / 0
                        Q : {(2,3),(7,6),(5,1),(4,8)}
                            no same edge
    (+1 -2 -3)(+4)  --> P1: {(2,4),(3,6),(5,1),(8,7)}  1-2-4-3-6-5-1, 7-8-7  1──>2───4<──3───6<──5──1, 7──>8──7  1 / 1
                        Q : {(2,3),(7,6),(5,1),(4,8)}
                                         ----- same edge
    (+1 -2 -4 -3)   --> P2: {(7,6),(5,1),(2,4),(3,8)}  1-2-4-3-8-7-6-5-1     1──>2───4<──3───8<──7───6<──5───1   1 / 2
                        Q : {(7,6),(5,1),(2,3),(4,8)}
                             ----- ----- same edges
    (+1 +2 -4 -3)   --> P3: {(7,6),(5,1),(2,3),(4,8)}  1-2-3-4-8-7-6-5-1     1──>2───3──>4───8<──7───6<──5───1   0 / 4
                        Q : {(7,6),(5,1),(2,3),(4,8)}
                             ----- ----- ----- ----- same edges
                             TODO: How can I implement 'BreakpointGraph(P,Q)' function?

Plan 1.
  * steps:
    a. randomly choose a blue-colored edge - (i,j)
    b. select two red-colored edges at the both ends of the blue edge (i-ip, j-jp)
    c. 2-BreakOnGenomeGraph(red, i,ip,j,jp)
       result: +1 trivial cycle
    d. convert modified set of red-colored edges to genome (GraphToGenome(P))
    e. add a new modifed genome to result collection

═════════════════════════════════════════════════

References:
- 582670 Algorithms for Bioinformatics
  https://www.cs.helsinki.fi/u/tpkarkka/teach/15-16/AfB/AfB_lecture5_20151001.pdf
- Best way to find the intersection of multiple sets?
  https://stackoverflow.com/questions/2541752/best-way-to-find-the-intersection-of-multiple-sets
- random.choice from set? python
  https://stackoverflow.com/questions/15837729/random-choice-from-set-python
- Python | Check if one list is subset of other
  https://www.geeksforgeeks.org/python-check-if-one-list-is-subset-of-other/
- How to clone or copy a set in Python?
  https://stackoverflow.com/questions/23200969/how-to-clone-or-copy-a-set-in-python
"""
#!/usr/bin/env python
import time
import random


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
    returns a chromosome from a cycle
    I don't think this is not simple function.
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])
        [1,-2,-3,4]
    >>> cycle_to_chromosome([2,4,3,6,5,1,2])
        [1,-2,-3]
    >>> cycle_to_chromosome([2,3,4,8,7,6,5,1,2])
        [2,-4,-3,1]
    """
    # if cycle[0] is cycle[-1], then remove either cycel[0] or cycle[-1]
    if cycle[0] == cycle[-1]:
        cycle = cycle[:-1]
    # In order to change the cycle to chromosome, it is necessary to first look at the cycle in detail.
    while not (cycle[0]%2 == 1 and cycle[1]%2 == 0 and (cycle[0]-1 == cycle[1] or cycle[0]+1 == cycle[1])):
        cycle = cycle[1:] + [cycle[0]]
    chromosome = []
    # as in textbook algorithm
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


def two_break_genome_graph(colored, i1, i2, j1, j2):
    """
    ({(int,int)},int,int,int,int) -> {(int,int)}
    returns a modifed colored edges
    As data structure of 'colored edges', I think the set is most appropriate.
    Blue-colored edges from genome Q never change.
    Red-colored edges from genome P will be modifed.
    The purpose of the 2-BreakOnGenomeGraph function is to make the genome P the same as the genome Q.
    >>> two_break_genome_graph({(2,4),(3,8),(7,5),(6,1)},1,6,3,8)
        {(2,4),(3,1),(6,8),(7,5)}
    """
    if (i1,i2) in colored:
        colored.remove((i1,i2))
        colored.add((i1,j1))
    else:
        colored.remove((i2,i1))
        colored.add((j1,i1))
    if (j1,j2) in colored:
        colored.remove((j1,j2))
        colored.add((i2,j2))
    else:
        colored.remove((j2,j1))
        colored.add((j2,i2))
    return colored


def colored_to_genome(colored):
    """
    {(a,a)} -> [[a]]
    returns a genome from colored edges
    This function modifies original set of colored edges. Be careful.
    >>> colored_to_genome({(2,4),(1,5),(6,7),(8,3)})
        [[2,-1,3,4]]
    """
    cycles = []
    cycle = []
    while colored:
        # if cycle is empty, start cycle
        if cycle == []:
            # pop an edge from a set of colored edges, init cycle
            cur = colored.pop()
            cycle += [cur[0], cur[1], cur[1]-1 if cur[1] % 2 == 0 else cur[1]+1]
        else:
            lst = cycle[-1]
            # find the next edge, and add nodes
            to_be_remove = ()
            for edge in colored:
                if edge[0] == lst:
                    cycle += [edge[1], edge[1]-1 if edge[1] % 2 == 0 else edge[1]+1]
                    to_be_removed = edge
                    break
                elif edge[1] == lst:
                    cycle += [edge[0], edge[0]-1 if edge[0] % 2 == 0 else edge[0]+1]
                    to_be_removed = edge
                    break
            colored.remove(to_be_removed)
        # if cycle closed, add cycle to cycles, and empty cycle
        if cycle[0] == cycle[-1]:
            cycles.append(cycle_to_chromosome(cycle))
            cycle = []
    return cycles


def two_break_on_genome(genome, io, ip, jo, jp):
    """
    ({(int,int)},int,int,int,int) -> [[int]]
    """
    colored = genome_to_colored(genome)
    twobreaks = two_break_genome_graph(colored,io,ip,jo,jp)
    genome = graph_to_genome(twobreaks)
    return genome


def breakpoint_graph(red, blue):
    """
    ({(a,a)},{(a,a)}) -> {(a,a)}
    returns a 'BreakPointGraph(red,blue)' from two sets of red and blue colored edges
    This function is not a simple union of two sets.
    >>> red  = {(2,3),(4,8),(5,1),(7,6)}
    >>> blue = {(2,4),(3,6),(5,1),(8,7)}   # <-- (5,1) overlaps in two sets.
    >>> red.union(blue)
        {(2,3),(4,8),(5,1),(7,6),(2,4),(3,6),(8,7)}  # one of two (5,1) was removed
    >>> breakpoint_graph(red,blue)
        {(2,3),(4,8),(5,1),(7,6),(2,4),(3,6),(5,1),(8,7)}  # <-- one of two (5,1) was substituted with (1,5)
    """
    common = set.intersection(red,blue)
    bpg = red.union(blue)
    for edge in common:
        bpg.add((edge[1],edge[0]))
    return bpg


def next_edge(cur, edges):
    """
    ((a,a),{(a,a)}) -> (a,a)
    returns the next edge following current edge, (cur,_) or (_,cur)
    >>> next_edge((1,2),{(2,1),(3,4)})
        (2,1)
    """
    next_edge = ()
    for edge in edges:
        if not cur:
            return next_edge
        if cur[0] in edge or cur[1] in edge:
            next_edge = edge
    return next_edge


def get_cycles_wrong(edges):
    """
    {(a,a)} -> [{a}]
    returns a list of cycles from a list of colored edges
    (don't care the order of elems in tuple edge)
    >>> p, q = [{(2,4),(3,6),(5,7),(8,1)}, {(2,3),(7,6),(5,1),(4,8)}]
    >>> get_cycles(p.union(q))
        [{(2,4),(4,8),(8,1),(5,1),(5,7),(7,6),(3,6),(2,3)}]
    >>> p, q = [{(2,4),(3,6),(8,7),(5,1)}, {(2,3),(7,6),(5,1),(4,8)}]
    >>> get_cycles(p.union(q))                   # notice the effect of union
        [{(2,4),(4,8),(8,7),(7,6),(3,6),(2,3)}]  # <-- trivial cycle removed automatically!
    >>> p, q = [{(2,3),(4,8),(7,6),(5,1)}, {(2,3),(7,6),(5,1),(4,8)}]  # all trivial cycles
    >>> get_cycles(p.union(q))
        []
    """
    cycles = []
    while edges:
        start = edges.pop()
        cycle = {start}
        cur = next_edge(start, edges)
        if not cur:
            break
        edges.remove(cur)
        while cur:
            cycle.add(cur)
            cur = next_edge(cur, edges)
            if not cur:
                break
            edges.remove(cur)
        cycles.append(cycle)
    return cycles


def nontrivial_cycles(edges):
    """
    {(a,a)} -> [{a}]
    returns a list of cycles from a list of colored edges
    (don't care the order of elems in tuple edge)
    >>> nontrivial_cycles({(1,9),(2,6),(3,2),(4,5),(5,12),(6,7),(7,3),(8,4),(9,8),(10,11),(11,10),(12,1)})
        [{(8,4),(9,8),(1,9),(12,1),(5,12),(4,5)},{(7,3),(3,2),(2,6),(6,7)}]
    >>> nontrivial_cycles({(2,4),(8,7),(5,1),(2,3),(7,6),(4,8),(3,6)})
    """
    cycles = []
    while edges:
        start = edges.pop()
        cycle = {start}
        cur = next_edge(start, edges)
        if not cur:
            break
        edges.remove(cur)
        while cur:
            cycle.add(cur)
            cur = next_edge(cur, edges)
            if not cur:
                break
            edges.remove(cur)
        if len(cycle) > 2:  # remove trivial cycles
            cycles.append(cycle)
    return cycles


def shortest_rearrangement_scenario(genome_p, genome_q):
    """
    ([[int]],[[int]]) -> [[int]]
    """
    result   = [genome_p]
    red    = genome_to_colored(genome_p)
    blue   = genome_to_colored(genome_q)
    cycles   = nontrivial_cycles(red.union(blue))
    while cycles:
        # select a cycle from cycles
        cycle = cycles[0]
        # select a blue edge
        i, j = random.choice(tuple(cycle.intersection(blue)))    # list(cycle - red)[0]
        # find two red-colored edges at the both ends of above blue-colored edge
        ip, jp = (-1, -1)
        for edge in red:
            if edge[0] == i or edge[1] == i:
                if edge[0] == i:
                    ip = edge[1]
                else:
                    ip = edge[0]
                break
        for edge in red:
            if edge[0] == j or edge[1] == j:
                if edge[0] == j:
                    jp = edge[1]
                else:
                    jp = edge[0]
                break
        # 2-BreakOnGenomeGraph
        red = two_break_genome_graph(red,i,ip,j,jp)
        # add modified genome P (P') to the result list
        genome = colored_to_genome(red.copy())
        result.append(genome)
        # update cycles
        cycles = nontrivial_cycles(red.union(blue))
    return result


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


def pretty_print_genome(genome):
    for chromosome in genome:
        str_chromosome = []
        for gene in chromosome:
            if gene >= 0:
                str_chromosome.append('+' + str(gene))
            else:
                str_chromosome.append(str(gene))
        print('(' + ' '.join(str_chromosome) + ')', end='')
    print()


def main():
    f = open('/home/wsl/rosalind/data/ba06d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    genome_p = read_genome(lines[0])
    genome_q = read_genome(lines[1])
    f.close()


    start_time = time.time()
    answer = shortest_rearrangement_scenario(genome_p, genome_q)
    answer
    for genome in answer:
        pretty_print_genome(genome)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
