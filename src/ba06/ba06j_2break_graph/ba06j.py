"""
Rosalind: BA6J
Implement 2-BreakOnGenomeGraph

We will use 2-Break(1, 6, 3, 8) two denote the 2-break replacing colored edges (1, 6) and (3, 8)
in a genome graph with two new colored edges (1, 3) and (6, 8).
Note that the order of the nodes in this function matter,
since the operation 2-Break(1, 6, 8, 3) would represent a different 2-break
that replaces (1, 6) and (3, 8) with (1, 8) and (6, 3).

The following pseudocode describes how 2-Break(i, i', j, j') transforms a genome graph.

╔═════════════════════════════════════════════════════════════════╗
║ 2-BreakOnGenomeGraph(GenomeGraph, i, i', j, j')                 ║
║      remove colored edges (i, i') and (j, j') from GenomeGraph  ║
║      add colored edges (i, j) and (i', j') to GenomeGraph       ║
║      return GenomeGraph                                         ║
╚═════════════════════════════════════════════════════════════════╝

2-Break On Genome Graph Problem
Solve the 2-Break On Genome Graph Problem.

Given: The colored edges of a genome graph GenomeGraph, followed by indices i, i', j, and j'.

Return: The colored edges of the genome graph resulting from applying the 2-break operation.

Sample Dataset
(2, 4), (3, 8), (7, 5), (6, 1)
1, 6, 3, 8

Sample Output
(2, 4), (3, 1), (7, 5), (6, 8)

═════════════════════════════════════════════════

    [ Where am I? ]

    * GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> Chromosome to Cycle (BA6F)
      -> Cycle to Chromosome (BA6G)
      -> Genome to ColoredEdges (BA6H)
      -> PREV: ColoredEdges to Genome (BA6I)
      -> HERE: 2-Break on ColoredEdges (BA6J)
      -> NEXT: 2-Break on Genome (BA6K)

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
  * just one collection (list of tuple)
    2-breaks   original    result
    -----------------------------
    (1,6,3,8)  (1,6)   --> (1,3)
               (3,8)   --> (6,8)
    (1,6,8,3)  (1,6)   --> (1,8)
               (3,8)   --> (6,3)

    where colored edge breaks: (1,6) (3,8)  same in both
                               (1,6,3,8)    and    (1,6,8,3)
    new colored edges        : (1,3) (6,8)         (1,8) (6,3)

  * a list of tuples
    breaks : [(i,i'),(j,j')]
    replace: (i,i')  --> (i,j)
             (i',i)  --> (j,i)
             (j,j')  --> (i',j')
             (j',j)  --> (j',i')

  * sample dataset
    (2,4), (3,8), (7,5), (6,1)      breaks: (1,6,3,8)
           -----         -----(rev)
    (2,4)  (6,8)  (7,5)  (3,1)      result 1
    (2,4)  (3,1)  (7,5)  (6,8)      result 2

    TODO: Are these results (1 and 2) the same? using BA6I
    colored_to_genome([(2,4),(6,8),(7,5),(3,1)]) == colored_to_genome([(2,4),(3,1),(7,5),(6,8)])
    teh same result: (+1 -2)(+3 -4)

  * regular expressions
    re.findall(r'\((.*?,.*?)\)', line)
    -- ------- - -------------
    a  b       c d
    (a) module 're' that stands for 'regular expression'
    (b) findall function
        https://docs.python.org/3/library/re.html#re.findall
    (c) symbol for starting regular expression 'r'
    (d) \(  : '(' symbol
        .   : any character except newline
        *   : 0 or more repetitions of the preceding RE
        ?   : match 0 or 1 repetitions of the preceding RE

Plan 1.
  * breaks = [(1,6), (3,8)]
  * remove all breaks --> add all new edges
    TODO: Does the order matter? No. Check the results of function in BA6I.

Plan 2.
  * breaks = [(1,6), (3,8)]
  * remove and add --> repeat

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time
import re


def two_break_genome_graph(colored, i1, i2, j1, j2):
    """
    ({(int,int)},int,int,int,int) -> {(int,int)}
    returns a modifed colored edges
    As data structure of 'colored edges', I think the set is most appropriate.
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


def line_to_graph(line):
    """
    str -> [(int,int)]
    returns a list of integer tuples (graph)
    >>> line_to_graph("(2,4),(3,8),(7,5),(6,1)")
        [(2, 4), (3, 8), (7, 5), (6, 1)]
    """
    graph = []
    ls_str = re.findall(r'\((.*?,.*?)\)', line)
    graph = set()
    for line in ls_str:
        fst, snd = [int(elem) for elem in line.split(', ')]
        graph.add((fst,snd))
    return graph


def main():
    f = open('/home/wsl/rosalind/data/ba06j.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    graph = line_to_graph(lines[0])
    i1, i2, j1, j2 = [int(elem) for elem in lines[1].split(', ')]
    f.close()

    start_time = time.time()
    answer = two_break_genome_graph(graph, i1, i2, j1, j2)
    print(', '.join(map(str, answer)))

    print(two_break_genome_graph({(2,4),(3,6),(5,7),(8,1)}, 2,4,3,6))
    print(two_break_genome_graph({(2,4),(3,6),(5,7),(8,1)}, 4,2,6,3))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()