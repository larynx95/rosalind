"""
Rosalind: BA6G
Implement CycleToChromosome

The process described in "Implement ChromosomeToCycle" is in fact invertible,
as described by the following pseudocode.

╔══════════════════════════════════════════════╗
║ CycleToChromosome(Nodes)                     ║
║      for j <- 1 to |Nodes| / 2               ║
║           if Node2j-1 < Node2j               ║
║                Chromosomej <- Node2j / 2     ║
║           else                               ║
║                Chromosomej <- -Node2j-1 / 2  ║
║      return Chromosome                       ║
╚══════════════════════════════════════════════╝

Cycle To Chromosome Problem
Solve the Cycle to Chromosome Problem.

Given:
A sequence Nodes of integers between 1 and 2n.

Return:
The chromosome Chromosome containing n synteny blocks resulting
from applying CycleToChromosome to Nodes.

Sample Dataset
(1 2 4 3 6 5 7 8)

Sample Output
(+1 -2 -3 +4)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> PREV: Chromosome to Cycle (BA6F)
      -> HERE: Cycle to Chromosome (BA6G)
      -> NEXT: Genome to ColoredEdge (BA6H)

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
  * sample dateset
    (1 2 4 3 6 5 7 8)
     --- --- --- ---
     ->  <-  <-  ->
     1+  -2  -3  +4

Plan 1.
  * algorithm
    ╔══════════════════════════════════════════════╗
    ║ CycleToChromosome(Nodes)                     ║
    ║      for j <- 1 to |Nodes| / 2               ║
    ║           if Node2j-1 < Node2j               ║
    ║                Chromosomej <- Node2j / 2     ║
    ║           else                               ║
    ║                Chromosomej <- -Node2j-1 / 2  ║
    ║      return Chromosome                       ║
    ╚══════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time


def cycle_to_chromosome(cycle):
    """
    [int] -> [int]
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])  # simple case
        [1,-2,-3,4]
    >>> cycle_to_chromosome([2,4,3,6,5,1,2])  # complicated case
        [1,-2,-3]
    >>> cycle_to_chromosome([2,3,4,8,7,6,5,1,2])  # more complicated case
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


def cycle_to_chromosome2(nodes):
    """
    [int] -> [str]
    returns a list of strings (same result as above)
    """
    chromosome = []
    for i in range(0, len(nodes), 2):
        p1, p2 = nodes[i], nodes[i+1]
        if p1 < p2:
            chromosome.append('+' + str(p2 // 2))
        else:
            chromosome.append('-' + str(p1 // 2))
    return chromosome


def cycle_to_chromosome3(nodes):
    """
    [int] -> [int]
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])
        [1,-2,-3,4]
    """
    Chromosome = []
    for i in range(0, len(nodes), 2):
        if nodes[i] < nodes[i + 1]:
            Chromosome.append(nodes[i + 1] // 2)
        else:
            Chromosome.append(-nodes[i] // 2)
    return Chromosome


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
    f = open('/home/wsl/rosalind/data/ba06g.txt', 'r')
    line = f.readline().strip('()')
    nodes = [int(elem) for elem in line.split(' ')]
    f.close()

    start_time = time.time()
    answer = cycle_to_chromosome(nodes)
    print(answer)
    pretty_print_genome([answer])
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
