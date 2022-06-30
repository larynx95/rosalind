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

References:
-
"""

def chromosome_to_cycle(chromosome):
    n = len(chromosome)
    nodes = [0] * 2 * n
    for j in range(n):
        i = int(chromosome[j])
        if i > 0:
            nodes[2 * j] = 2 * i - 1
            nodes[2 * j + 1] = 2 * i
        else:
            nodes[2 * j] = -2 * i
            nodes[2 * j + 1] = -2 * i - 1
    return nodes


def colored_edges(p):
    edges = []
    for chromosome in p:
        nodes = chromosome_to_cycle(chromosome)
        nodes.append(nodes[0])
        for i in range(len(chromosome)):
            edge = (nodes[2*i+1], nodes[2*i+2])
            edges.append(edge)
    return edges


def main():
    p = input().split(')(')
    p[0], p[-1] = p[0][1:], p[-1][:-1]
    p = list(map(lambda x: x.split(), p))
    c = list(map(str, colored_edges(p)))
    print(', '.join(c))


if __name__ == '__main__':
    main()