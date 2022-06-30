"""
Rosalind: BA6H
Implement ColoredEdges

The following algorithm constructs ColoredEdges(P) for a genome P.
In this pseudocode,
we will assume that an n-element array (a1, . . . , an) has an invisible (n + 1)-th element
that is equal to its first element, i.e., an+1 = a1.

╔═════════════════════════════════════════════════════════════╗
║ ColoredEdges(P)                                             ║
║      Edges <- an empty set                                  ║
║      for each chromosome Chromosome in P                    ║
║           Nodes <- ChromosomeToCycle(Chromosome)            ║
║           for j <- 1 to |Chromosome|                        ║
║                add the edge (Nodes2j, Nodes2j +1) to Edges  ║
║      return Edges                                           ║
╚═════════════════════════════════════════════════════════════╝

Colored Edges Problem
Find the Colored Edges in a genome.

Given: A genome P.

Return: The collection of colored edges in the genome graph of P in the form (x, y).

Sample Dataset
(+1 -2 -3)(+4 +5 -6)

Sample Output
(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> PREV: Cycle to Chromosome (BA6G)
      -> HERE: Genome to ColoredEdge (BA6H)
      -> NEXT: ColoredEdge to Genome (BA6I)

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
  * TODO: Does the order of colored edges matter?
          No, the order of colored edges doesn't matter.

  * chromosome to nodes - (+1 -2 -3)(+4 +5 -6)
    (1) the first chromosome (+1 -2 -3)
                            +a      -b       -c
    linear chromosome    : <───*───*───>*───*───>
    -> permutation       : (+1      -2       -3)
    -> directed edges    : <-1       2->      3->
    -> undirected edges  : (1t-1h   2h-2t    3h-3t)
    -> nodes (integers)  :  1  2    4  3     6  5    1
                               ------  -------  ------
                               (2, 4)  (3, 6)   (5, 1)

    (2) the second chromosome (+4 +5 -6)
                            +a      +b       -c
    linear chromosome    : <───*───*───>*───*───>
    -> permutation       : (+4     +5     -6)
    -> directed edges    : <-4    <-5      6->
    -> undirected edges  : (4t-4h  5t-5h  6h-6t)
    -> nodes (integers)  :  4  5   6  7   9  8     incorrect! overlapping integers

                            +a      +b       -c
    linear chromosome    : <───*───*───>*───*───>
    -> permutation       : (+8     +9     -10)
    -> directed edges    : <-8    <-9      10->
    -> undirected edges  : (8t-8h  9t-9h  10h-10t)
    -> nodes (integers)  :  8  9  10 11   13  12   incorrect!

                            +a      +b       -c
    linear chromosome    : <───*───*───>*───*───>
    -> permutation       : (+7     +9     -11)
    -> directed edges    : <-7    <-9      11->
    -> undirected edges  : (7t-7h  9t-9h  11h-11t)
    -> nodes (integers)  :  7  8   9  10  12  11    7   correct
                               -----  ------  -------
                               (8, 9) (10, 12) (11, 7)

  * nodes to colored edges
    n    0       1       2        0       1         2
    2n   0       2       4        0       2         4
    2n+1    1       3       5        1        3         5
    idx  0  1    2  3    4  5     0  1    2   3     4   5
         1  2    4  3    6  5     7  8    9  10    12  11

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time


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


def genome_to_colored1(genome):
    """
    [[int]] -> {(int,int)}
    returns a list of list of integers
    >>> genome_to_colored1([[1,-2,-3],[4,5,-6]])
        {(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)}
    """
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


def genome_to_colored2(P):
    edges = list()
    for chromosome in P:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(1, len(nodes), 2):
            if j != len(nodes) - 1:
                edges.append((nodes[j], nodes[j + 1]))
            else:
                edges.append((nodes[j], nodes[0]))
    return edges


def genome_to_colored3(genome):
    """
    [[int]] -> [(int,int)]
    returns a list of list of integers
    >>> genome_to_colored3([[1,-2,-3],[4,5,-6]])
        [(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)]
    """
    edges = []
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.append((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


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
    f = open('/home/wsl/rosalind/data/ba06h.txt', 'r')
    strings = f.readline().strip('\n')[1:-1].split(')(')
    chromosomes = []
    for string in strings:
        chromosomes.append([int(elem) for elem in string.split(' ')])
    f.close()

    start_time = time.time()
    # This answer is correct.
    answer = [str(tup) for tup in genome_to_colored1(chromosomes)]
    print(', '.join(answer))

    # But what if genes are not in order? These are incorrect.
    print(genome_to_colored3([[+1,-2,-3,+4]]))  # [(2,4),(3,6),(5,7),(8,1)]
    print(genome_to_colored3([[+1,+2,-4,-3]]))  # [(2,3),(4,8),(7,6),(5,1)]
    print(genome_to_colored3([[+4]]))           # [(8,7)]
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
