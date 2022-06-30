"""
Rosalind: BA6I
Implement GraphToGenome

The colored edges in the breakpoint graph of P and Q are given by ColoredEdges(P) together with ColoredEdges(Q).
Note that some edges in these two sets may connect the same two nodes, which results in trivial cycles.

We will find it helpful to implement a function converting a genome graph back into a genome.
╔═══════════════════════════════════════════════════╗
║ GraphToGenome(GenomeGraph)                        ║
║      P <- an empty set of chromosomes             ║
║      for each cycle Nodes in GenomeGraph          ║
║           Chromosome <- CycleToChromosome(Nodes)  ║
║           add Chromosome to P                     ║
║      return P                                     ║
╚═══════════════════════════════════════════════════╝

Graph To Genome Problem
Solve the Graph To Genome Problem.

Given: The colored edges of a genome graph.

Return: A genome corresponding to the genome graph.

Sample Dataset
(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)

Sample Output
(+1 -2 -3)(-4 +5 -6)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> Cycle to Chromosome (BA6G)
      -> PREV: Genome to ColoredEdges (BA6H)
      -> HERE: ColoredEdges to Genome (BA6I)
      -> NEXT: 2-Break on ColoredEdges (BA6J)

    [ Summary of Helper Functions ]

    * Terminology         :: Type
      genome              :: [[int]]
      chromosome          :: [int]
      cycle/nodes         :: [int]
      colored/black edges :: {(int,int)}
      graph
        adjacency list    :: {(a,b)}, [(a,b)], [(a,[b])]
        dictionary        :: {a:b}, {a:[b]}

    * Functions
      BA6F - ChromosomeToCycle
      BA6G - CycleToChromosome
      BA6H - ColoredEdges (GenomeToColoredEdges)
      BA6I - GraphToGenome (ColoredEdgesToGenome)
      BA6J - 2-BreakOnGenomeGraph
      BA6K - 2-BreakOnGenone

      ChromosomeToCycle  -> ColoredEdges ─┬─┐
                    2-BreakOnGenomeGraph ─┘ │
      CycleToChromosome -> GraphToGenome ───┴─── 2-BreakOnGenome

Info.
  * TODO: Does the order matter?
    case 1. (2,4),(3,6),(5,1),(8,9),(10,12),(11,7)
           1->2-4->3-6->5-1  7->8-9->10-12->11-7
           +1   -2   -3      +4   +5    -6

    case 2. (2,4),(3,6),(11,7),(5,1),(8,9),(10,12)
           1->2-4->3-6->5-1   11-7->8-9->10-12->11
           +1   -2   -3          +4   +5    -6

    case 3. (2,4),(8,9),(3,6),(10,12),(11,7),(5,1)
           1->2-4->3-6->5-1   10-12->11-7->8-9->10
           +1   -2   -3          -6     +4   +5

    No, the order of colored edges doesn't matter.
    So the 'GraphToGenome' function must return the same result
    regardless of whether or not the arguments are aligned.

  * sample dataset
    colored edges: (2,4),(3,6),(5,1),(7,9),(10,12),(11,8)
              1  2  3  4  5  6  7  8  9 10 11 12
              ----  ----  ----  ----  ---- -----    directed, black edges, synteny block
    ( 2, 4): (1  2  4  3
    ( 3, 6): (1  2  4  3  6  5
    ( 5, 1): (1  2  4  3  6  5)
    ( 7, 9): (1  2  4  3  6  5)(8  7  9 10
    (10,12): (1  2  4  3  6  5)(8  7  9 10 12 11
    (11, 8): (1  2  4  3  6  5)(8  7  9 10 12 11)
              ----  ----  ----  ----  ---- -----
             ( +1    -2    -3 )( -4    +5   -6 )

Plan 1.
  * one step
  * genome: (+1 -2 -3) (-4 +5 -6)
    cycle                         colored edges (graph)
    []                            (10,12) (5,1) (2,4) (3,6) (7,9) (11,8)
    [9,10,12,11]                          (5,1) (2,4) (3,6) (7,9) (11,8)
    [9,10,12,11,8,7]                      (5,1) (2,4) (3,6) (7,9)
    [9,10,12,11,8,7,*9*]                  (5,1) (2,4) (3,6)
     (+5   -6    -4)  <-- first chromosome
    []                                    (5,1) (2,4) (3,6)
    [6,5,1,2]                                   (2,4) (3,6)
    [6,5,1,2,4,3]                                     (3,6)
    [6,5,1,2,4,3,*6*]
    (-3  +1  -2)  <-- second chromosome

Plan 2.
  * two steps
    - merge colored and black edges (Edges create cycles.)
    - get cycles from the set (and cycle to chromosome)

═════════════════════════════════════════════════

References:
- Python Regular expression matching tuple pairs
  https://stackoverflow.com/questions/48030275/python-regular-expression-matching-tuple-pairs
- Set changed size during iteration
  https://stackoverflow.com/questions/24610784/set-changed-size-during-iteration
- List of tuples to dictionary [duplicate]
  https://stackoverflow.com/questions/6522446/list-of-tuples-to-dictionary
- Access an arbitrary element in a dictionary in Python
  https://stackoverflow.com/questions/3097866/access-an-arbitrary-element-in-a-dictionary-in-python
- Accessing dict_keys element by index in Python3
  https://stackoverflow.com/questions/18552001/accessing-dict-keys-element-by-index-in-python3
"""
#!/usr/bin/env python
import time
import re


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


def colored_to_black_wrong(colored):
    """
    {(a,a)} -> {(a,a)}
    return a set of black edges from colored edges
    Each black edge is directed edge. So this function is wrong
    >>> colored_to_black_wrong([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        {(9,10),(1,2),(6,5),(4,3),(8,7),(12,11)}
    """
    black = set()
    for edge in colored:
        if edge[0] % 2 == 0:
            black.add((edge[0]-1,edge[0]))
        else:
            black.add((edge[0]+1,edge[0]))
        if edge[1] % 2 == 0:
            black.add((edge[1],edge[1]-1))
        else:
            black.add((edge[1],edge[1]+1))
    return black


def colored_to_genome_wrong(colored):
    """
    {(int,int)} -> [[int]]
    returns a genome from colored edges (BA6I)
    >>> colored_to_genome_wrong({(10,12),(5,1),(2,4),(3,6),(7,9),(11,8)})
        [[1,-2,-3],[-4,5,-6]]
    sometimes works, sometimes not TODO: Fix this.
    """
    genome = []
    cycle = []
    cur = ()
    while colored:
        # initialize a cycle
        if cycle == []:
            cur = colored.pop()
            cycle = [cur[0] - 1 if cur[0] % 2 == 0 else cur[0] + 1, cur[0],
                     cur[1], cur[1] - 1 if cur[1] % 2 == 0 else cur[1] + 1]
        else:
            # find next edge (to avoid 'Set changed size during iteration' error)
            for edge in colored:
                if cycle[-1] == edge[0]:
                    cur = edge
            # condition for end of cycle
            if cur[1] == cycle[0]:
                genome.append(cycle_to_chromosome(cycle))
                cycle = []
            else:
                cycle += [cur[1]] + [cur[1] -1 if cur[1] % 2 == 0 else cur[1] + 1]
            colored.remove(cur)
    return genome


def colored_to_genome_wrong2(colored):
    """
    {(a,a)} -> [[a]]
    returns a genome from colored edges
    >>> colored_to_genome_wrong2({(2, 4), (1, 5), (6, 7), (8, 3)})
        error
    """
    # get a set of black edges
    black = set()
    for edge in colored:
        if edge[0] % 2 == 0:
            black.add((edge[0]-1,edge[0]))
        else:
            black.add((edge[0]+1,edge[0]))
        if edge[1] % 2 == 0:
            black.add((edge[1],edge[1]-1))
        else:
            black.add((edge[1],edge[1]+1))
    # merge colored and black edges into a new dictionary
    dic = dict(colored.union(black))
    genome = []
    cycle = []
    cur = -1
    while dic:
        if cycle == []:
            cur = list(dic)[0]
            val = dic[cur]
            cycle += [cur, val]
            del dic[cur]
            cur = val
        elif cur not in dic:
            genome.append(cycle_to_chromosome(cycle))
            cycle = []
        else:
            val = dic[cur]
            cycle += [val]
            del dic[cur]
            cur = val
        if not dic:
            genome.append(cycle_to_chromosome(cycle))
    return genome


def colored_to_genome(colored):
    """
    {(a,a)} -> [[a]]
    returns a genome from colored edges
    This function modifies original set of colored edges. Be careful.
    >>> colored_to_genome_wrong3({(2,4),(1,5),(6,7),(8,3)})
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


def pretty_print_genome(genome):
    for chromosome in genome:
        str_chromosome = []
        for gene in chromosome:
            if gene >= 0:
                str_chromosome.append('+' + str(gene))
            else:
                str_chromosome.append(str(gene))
        print('(' + ' '.join(str_chromosome) + ')', end=' ')
    print()


def main():
    f = open('/home/wsl/rosalind/data/ba06i.txt', 'r')
    line = f.readline().strip('\n')
    temp = re.findall(r'\((.*?,.*?)\)', line)
    colored = set()
    for frag in temp:
        fst, snd = [int(elem) for elem in frag.split(',')]
        colored.add((fst,snd))
    f.close()

    start_time = time.time()
    genome = colored_to_genome(colored)
    pretty_print_genome(genome)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
