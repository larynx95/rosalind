"""
Rosalind: BA6F
Implement ChromosomeToCycle

The following pseudocode bypasses the intermediate step of assigning "head" and "tail" nodes
in order to transform a single circular chromosome
Chromosome = (Chromosome1, . . . , Chromosomen)
into a cycle represented as a sequence of integers
Nodes = (Nodes1, . . . , Nodes2n).

╔════════════════════════════════════╗
║ ChromosomeToCycle(Chromosome)      ║
║      for j <- 1 to |Chromosome|    ║
║           i <- Chromosomej         ║
║           if i > 0                 ║
║                Node2j-1 <- 2i-1    ║
║                Node2j   <- 2i      ║
║           else                     ║
║                Node2j-1 <-  -2i    ║
║                Node2j <- -2i-1     ║
║      return Nodes                  ║
╚════════════════════════════════════╝

Chromosome To Cycle Problem

Solve the Chromosome To Cycle Problem.

Given:
A chromosome Chromosome containing n synteny blocks.

Return:
The sequence Nodes of integers between 1 and 2n
resulting from applying ChromosomeToCycle to Chromosome.

Sample Dataset
(+1 -2 -3 +4)

Sample Output
(1 2 4 3 6 5 7 8)

═════════════════════════════════════════════════

    [ Where am I? ]

    * PREV: GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> HERE: ChromosomeToCycle (BA6F)
      -> NEXT: CycleToChromosome (BA6G)

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
  * sample dataset
    (+1     -2     -3     +4)
     1t->1h 2t->2t 3h->3t 4h->4t
     1 -> 2 4 <- 3 6 <- 5 7 -> 8
    ( 1  2   4  3   6  5   7  8 )

Plan 1.
  * The way how chromosome is represented in this chapter is something I never thought of.
    That was the most difficult thing to me in this chapter.

  * representation of chromosome (fig. 6-24)
                            +a      -b       -c       +d
    linear chromosome    : <───*───*───>*───*───>*───*<───
    -> permutation       : (+1     -2     -3     +4)         signed integers
    -> directed edges    : <-1      2->    3->  <-4           circular graph
    -> undirected edges  : (1t-1h  2h-2t  3h-3t  4t-4h)      tail(t), head(h), t->h(+), h->t(-)
    -> nodes (integers)  :  1  2   4  3   6  5   7  8

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
    f = open('/home/wsl/rosalind/data/ba06f.txt', 'r')
    line = f.readline().strip('()')
    perm = [int(elem) for elem in line.split(' ')]
    f.close()

    start_time = time.time()
    answer = list(map(str, chromosome_to_cycle(perm)))
    print('(' + ' '.join(answer) + ')')
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()