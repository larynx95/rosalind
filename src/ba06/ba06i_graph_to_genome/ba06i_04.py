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
import re


def cycle_to_chromosome(nodes):
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


def graph_to_genome(graph):
    """
    [(int,int)] -> [[int]]
    >>> graph_to_genome([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        [[1,-2,-3],[-4,5,-6]]
    """
    P = []
    Cycles = []
    temp = []
    for i in range(len(graph)):
        if i == len(graph) - 1:
            temp += graph[i]
            Cycles.append(temp)
        elif graph[i][1] == graph[i + 1][0] + 1 or graph[i][1] == graph[i + 1][0] - 1:
            temp += graph[i]
        else:
            temp += graph[i]
            Cycles.append(temp)
            temp = []
    for Cycle in Cycles:
        Chromosome = cycle_to_chromosome([Cycle[-1]] + Cycle[:-1])
        P.append(Chromosome)
    return P


result1 = graph_to_genome([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
result2 = graph_to_genome([(7,9),(10,12),(2,4),(3,6),(5,1),(11,8)])
result3 = graph_to_genome([(5,1),(2,4),(3,6),(7,9),(11,8),(10,12)])
print(result1 == result2 == result3)  # False
print(graph_to_genome([(2,4),(3,8),(7,5),(6,1)]))  # [['+1','-2','-4'],['+3']]


def main():
    f = open('/home/wsl/rosalind/data/ba06i.txt', 'r')
    line = f.readline().strip('\n')
    temp = re.findall(r'\((.*?,.*?)\)', line)
    nodes = []
    for frag in temp:
        fst, snd = [int(elem) for elem in frag.split(',')]
        nodes.append((fst,snd))
    f.close()

    result = graph_to_genome(nodes)
    for j in range(len(result)):
        result[j] = '(' + ' '.join(('+' if i > 0 else '') + str(i) for i in result[j]) + ')'
    print(''.join(result))


if __name__ == '__main__':
    main()
