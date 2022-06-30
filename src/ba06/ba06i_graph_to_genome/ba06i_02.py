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

References:
-
"""
#!/usr/bin/env python
def cycle_to_chromosome(nodes):
    """
    [int] -> [str]
    >>> cycle_to_chromosome([1,2,4,3,6,5])
        ['+1','-2','-3']
    """
    chromosome = []
    for i in range(0, len(nodes), 2):
        p1, p2 = nodes[i], nodes[i+1]
        if p1 < p2:
            chromosome.append('+' + str(p2 // 2))
        else:
            chromosome.append('-' + str(p1 // 2))
    return chromosome


def graph_to_genome(genome_graph):
    """
    [(int,int)] -> str
    >>> graph_to_genome([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        (+1 -2 -3)(-4 +5 -6)
    >>> graph_to_genome([(2, 4), (3, 8), (7, 5), (6, 1)])
        (+1 -2 -4 +3)
    """
    p = []
    nodes = []
    for edge in genome_graph:
        adj_left = edge[0] + 1 if edge[0] % 2 == 1 else edge[0] - 1
        nodes += [adj_left, edge[0]]
        if nodes[0] == edge[1]:
            p += '(' + ' '.join(cycle_to_chromosome(nodes)) + ')'
            nodes = []
    return p


def main():
    f = open('/home/wsl/rosalind/data/ba06i.txt', 'r')
    line = f.readline().strip('\n')
    f.close()
    # genome_graph = eval('[' + input() + ']')
    genome_graph = eval('[' + line + ']')
    print(graph_to_genome(genome_graph))


if __name__ == '__main__':
    main()
