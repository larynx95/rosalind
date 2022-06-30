"""
Rosalind: BA6K
Implement 2-BreakOnGenome

We can extend the pseudocode in "Implement 2-BreakOnGenomeGraph" to a 2-break defined on genome P.

╔══════════════════════════════════════════════════════════════════════╗
║ 2-BreakOnGenome(P, i, i', j, j')                                     ║
║      GenomeGraph <- BlackEdges(P) and ColoredEdges(P)                ║
║      GenomeGraph <- 2-BreakOnGenomeGraph(GenomeGraph, i, i', j, j')  ║
║      P <- GraphToGenome(GenomeGraph)                                 ║
║      return P                                                        ║
╚══════════════════════════════════════════════════════════════════════╝

Implement 2-BreakOnGenome
Solve the 2-Break On Genome Graph Problem.

Given: A genome P, followed by indices i, i', j, and j'.

Return: The genome P' resulting from applying the 2-break operation.

Sample Dataset
(+1 -2 -4 +3)
1, 6, 3, 8

Sample Output
(+2 -1) (-3 +4)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> Cycle to Chromosome (BA6G)
      -> Genome to ColoredEdges (BA6H)
      -> ColoredEdges to Genome (BA6I)
      -> PREV: 2-Break on ColoredEdges (BA6J)
      -> HERE: 2-Break on Genome (BA6K)
      ↓
    * NEXT: 2-Break distance (BA6C)

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
  * Sample Dataset
    genome          : (+1   -2   -4   +3)
    nodes           :   1  2 4  3 8  7 5  6
    colored         :     (2 4)(3 8)(7 5)(6 1)
    break           : (1 6)(3 8)
    after 2-break   :     (2 4)(3 1)(7 5)(6 8)
    genome          : (-2 +1) (+3 -4)

Plan 1.
  * steps: BA6H -> BA6J -> BA6I -> BA6K
    (1) genome        to colored edges         - BA6H
    (2) colored edges to 2-break colored edges - BA6J
    (3) colored edges to genome                - BA6I
    All the functions in previous exercises are correct.
    But the combined result from those functions is not correct. Why?
    Sometimes it's right and sometimes it's not. (TODO: Why?)

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time
import re


def chromosome_to_cycle(chromosome):
    """
    [int] -> [int]
    algorithm in textbook (BA6F)
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
    returns a set of list of integers (BA6H)
    >>> genome_to_colored([[1,-2,-3],[4,5,-6]])
        {(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)}
    >>> genome_to_colored([[+1,-2,-4,+3]])
        {(6,1),(2,4),(7,5),(3,8)}
    """
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


def genome_to_black(genome):
    """
    [[int]] -> {(int,int)}
    returns a set of black edges
    >>> genome_to_black([[+1,-2,-4,+3]])
        {(1,2),(4,3),(8,7),(5,6)}
    """
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j], nodes[2*j+1 % len(nodes)]))
    return edges


def two_break_genome_graph(colored, i1, i2, j1, j2):
    """
    ({(int,int)},int,int,int,int) -> {(int,int)}
    returns a modifed colored edges (BA6J)
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
    # colored edges to genome
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


def two_break_on_genome(genome, i1, i2, j1, j2):
    """
    ({(int,int)},int,int,int,int) -> [[int]]
    """
    colored = genome_to_colored(genome)
    twobreaks = two_break_genome_graph(colored,i1,i2,j1,j2)
    genome = colored_to_genome(twobreaks)
    return genome


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
    f = open('/home/wsl/rosalind/data/ba06k.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    chromosome = [int(elem) for elem in lines[0][1:-1].split(' ')]
    i1, i2, j1, j2 = [int(elem) for elem in lines[1].split(', ')]
    f.close()

    start_time = time.time()
    # multiple steps
    #colored = genome_to_colored([chromosome])
    #twobreaks = two_break_genome_graph(colored,i1,i2,j1,j2)
    #genome = colored_to_genome(twobreaks)
    #pretty_print_genome(genome)

    # one step
    genome = two_break_on_genome([chromosome],i1,i2,j1,j2)
    pretty_print_genome(genome)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
