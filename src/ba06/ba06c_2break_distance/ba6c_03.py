"""
Rosalind: BA6C
Compute the 2-Break Distance Between a Pair of Genomes

2-Break Distance Problem
Find the 2-break distance between two genomes.

Given: Two genomes with circular chromosomes on the same set of synteny blocks.

Return: The 2-break distance between these two genomes.

Sample Dataset
(+1 +2 +3 +4 +5 +6)
(+1 -3 -6 -5)(+2 -4)

Sample Output
3

═════════════════════════════════════════════════

References:
* solution by other person
- https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA6C.py
"""
#!/usr/bin/env python
import sys


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
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])
        [1,-2,-3,4]
    >>> cycle_to_chromosome([2,4,3,6,5,1,2])  # <-- prepare for extreme case
    """
    if cycle[0] == cycle[-1]:
        cycle = cycle[:-1]
    while cycle[0]+1 != cycle[1] and cycle[0]-1 != cycle[1]:
        cycle = cycle[1:] + [cycle[0]]
    chromosome = []
    for i in range(0, len(cycle), 2):
        if cycle[i] < cycle[i + 1]:
            chromosome.append(cycle[i + 1] // 2)
        else:
            chromosome.append(-cycle[i] // 2)
    return chromosome


def colored_edges(P):
    Edges = list()
    for chromosome in P:
        Nodes = chromosome_to_cycle(chromosome)
        for j in range(1, len(Nodes), 2):
            if j != len(Nodes) - 1:
                Edges.append([Nodes[j], Nodes[j + 1]])
            else:
                Edges.append([Nodes[j], Nodes[0]])
    return Edges


def find_next_edge(current, edges):
    if len(edges) == 0:
        return -1
    idx = 0
    while not (current[0] in edges[idx] or current[1] in edges[idx]):
        idx += 1
        if idx == len(edges):
            return -1
    return edges[idx]


def two_break_distance(P, Q):
    edgesP = colored_edges(P)
    edgesQ = colored_edges(Q)
    edges = edgesP + edgesQ
    blocks = set()
    for edge in edges:
        blocks.add(edge[0])
        blocks.add(edge[1])
    Cycles = []
    while len(edges) != 0:
        start = edges[0]
        edges.remove(edges[0])
        Cycle = [start]
        current = find_next_edge(start, edges)
        while current != -1:
            Cycle.append(current)
            edges.remove(current)
            current = find_next_edge(current, edges)
        Cycles.append(Cycle)
    return len(blocks) // 2 - len(Cycles)


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


if __name__ == "__main__":
    '''
    Given: Two genomes with circular chromosomes on the same set of synteny blocks.
    Return: The 2-break distance between these two genomes.
    '''
    f = open('/home/wsl/rosalind/data/ba06c.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    genome_p = read_genome(lines[0])
    genome_q = read_genome(lines[1])
    f.close()

    print(two_break_distance(genome_p, genome_q))
