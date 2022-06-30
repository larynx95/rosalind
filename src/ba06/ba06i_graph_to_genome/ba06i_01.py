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

collection of incorrect functions
"""
#!/usr/bin/env python
import time
import re


def chromosome_to_cycle(chromosome):
    """
    [int] -> [int]
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


def cycle_to_chromosome(nodes):
    """
    [int] -> [str]
    return a list of strings (signed integers)
    >>> cycle_to_chromosome([1,2,4,3,6,5,7,8])
        (+1 -2 -3 +4)
    """
    chromosome = []
    for j in range(0, len(nodes) // 2):
        if nodes[2*j] < nodes[2*j + 1]:
            chromosome.append('+' + str(nodes[2*j + 1] // 2))
        else:
            chromosome.append('-' + str(nodes[2*j] // 2))
    return chromosome


def genome_to_colored(chromosomes):
    """
    [[int]] -> {(int,int)}
    returns a list of list of integers
    >>> genome_to_colored([[1,-2,-3],[4,5,-6]])
        {(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)}
    """
    edges = set()
    for chromosome in chromosomes:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


def colored_to_black(list_colored_edges):
    """
    [[(int,int)]] -> [[int]]
    returns a list of black edges from colored edges
    useless function, ignore this function
    >>> colored_to_black([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        [(1,2),(4,3),(6,5),(8,7),(9,10),(12,11)]
    """
    nodes = set()
    for edge in list_colored_edges:
        if edge[0] % 2 == 0:
            nodes.add((edge[0] - 1, edge[0]))
        else:
            nodes.add((edge[0] + 1, edge[0]))
        if edge[1] % 2 == 0:
            nodes.add((edge[1], edge[1] - 1))
        else:
            nodes.add((edge[1], edge[1] + 1))
    return sorted(list(nodes))


def black_to_cycle(black_edges):
    """
    [(int,int)] -> [int]
    returns a list of integer list
    (Problem!) Black edges can't construct multiple chromosomes.
    >>> black_to_cycle([(1,2),(4,3),(6,5),(8,7),(9,10),(12,11)])
        [1,2,4,3,6,5,8,7,9,10,12,11]
    """
    nodes = []
    for tup in black_edges:
        nodes.append(tup[0])
        nodes.append(tup[1])
    return nodes


def colored_to_cycles(colored_edges):
    """
    [(int,int)] -> [[int]]
    returns a list of integer lists (lists of nodes, cycles)
    >>> colored_to_cycles([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        [[1,2,4,3,6,5],[8,7,9,10,12,11]]
    """
    colored_edges = sorted(colored_edges)  # <-- Is this right?
    cycles = []
    temp = []
    for edge in colored_edges:
        temp.append(edge[0])
        if edge[0] > edge[1]:
            temp = [edge[1]] + temp
            cycles.append(temp)
            temp = []
        else:
            temp.append(edge[1])
    return cycles


#################################################
# case 1
#################################################

def graph_to_genome(colored_edges):
    """
    [(int,int)] -> [str]
    implementation of 'GraphToGeneome' function
    (Caution) This function is not always correct! TODO: Fix this!
    >>> graph_to_genome([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        [['+1','-2','-3'],['-4','+5','-6']]
    >>> graph_to_genome([(2,4),(3,8),(7,5),(6,1)])  # (+1 -2 -4 +3)
        [['+1','-2','-4'],['+3']]    # <-- incorrect! TODO: Fix this!
    """
    colored_edges = sorted(colored_edges)  # <-- Is this right?
    cycles = []
    temp = []
    for edge in colored_edges:
        temp.append(edge[0])
        if edge[0] > edge[1]:
            temp = [edge[1]] + temp
            cycles.append(temp)
            temp = []
        else:
            temp.append(edge[1])
    chromosomes = []
    for cycle in cycles:
        chromosomes.append(cycle_to_chromosome(cycle))
    return chromosomes


result1 = graph_to_genome([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
result2 = graph_to_genome([(7,9),(10,12),(2,4),(3,6),(5,1),(11,8)])
result3 = graph_to_genome([(5,1),(2,4),(3,6),(7,9),(11,8),(10,12)])
print(result1 == result2 == result3)  # True
print(graph_to_genome([(2,4),(3,8),(7,5),(6,1)]))  # [['+1','-2','-4'],['+3']]


#################################################
# case 2
#################################################

def graph_to_genome2(genome_graph):
    """
    [(int,int)] -> str
    >>> graph_to_genome2([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        (+1 -2 -3)(-4 +5 -6)
    >>> graph_to_genome2([(2, 4), (3, 8), (7, 5), (6, 1)])
        (+1 -2 -4 +3)
    """
    p = ''
    nodes = []
    for edge in genome_graph:
        adj_left = edge[0] + 1 if edge[0] % 2 == 1 else edge[0] - 1
        nodes += [adj_left, edge[0]]
        if nodes[0] == edge[1]:
            p += '(' + ' '.join(cycle_to_chromosome(nodes)) + ')'
            nodes = []
    return p


result1 = graph_to_genome2([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
result2 = graph_to_genome2([(7,9),(10,12),(2,4),(3,6),(5,1),(11,8)])
result3 = graph_to_genome2([(5,1),(2,4),(3,6),(7,9),(10,12),(11,8)])
print(result1 == result2 == result3)  # False
print(graph_to_genome2([(2,4),(3,8),(7,5),(6,1)]))  # (+1 -2 -4 +3)

#################################################
# case 3
#################################################

def graph_to_genome3(genome_graph):
    """
    [(int,int)] -> str
    >>> graph_to_genome3([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
        [['+1','-2','-3'],['-4','+5','-6']]
    >>> graph_to_genome3([(2,4),(3,8),(7,5),(6,1)])
        [['+1','-2','-4','+3']]
    """
    genome = []
    nodes = []
    for edge in genome_graph:
        adj_left = edge[0] + 1 if edge[0] % 2 == 1 else edge[0] - 1
        nodes += [adj_left, edge[0]]
        if nodes[0] == edge[1]:
            genome.append(cycle_to_chromosome(nodes))
            nodes = []
    return genome


result1 = graph_to_genome3([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
result2 = graph_to_genome3([(7,9),(10,12),(2,4),(3,6),(5,1),(11,8)])
result3 = graph_to_genome3([(5,1),(2,4),(3,6),(7,9),(10,12),(11,8)])
print(result1 == result2 == result3)  # False
print(graph_to_genome3([(2,4),(3,8),(7,5),(6,1)]))  # [['+1', '-2', '-4', '+3']]


#################################################
# case 4
#################################################

def FindNextEdge(current, edges):
    if len(edges) == 0:
        return -1
    idx = 0
    while not (current[1] + 1 == edges[idx][0] or current[1] - 1 == edges[idx][0]):
        idx += 1
        if idx == len(edges):
            return -1
    return edges[idx]


def graph_to_genome4(GenomeGraph):
    """
    >>> graph_to_genome4([[2,4],[6,8],[7,5],[3,1]])
        [[-2,1],[-4,3]]
    """
    Q = []
    Cycles = []
    idx = 0
    while len(GenomeGraph) != 0:
        Cycle = []
        current = GenomeGraph[0]
        while current != -1:
            Cycle += current
            GenomeGraph.remove(current)
            current = FindNextEdge(current, GenomeGraph)
        Cycles.append(Cycle)
    for Cycle in Cycles:
        Cycle = Cycle[-3:] + Cycle[:-3]
        Chromosome = cycle_to_chromosome(Cycle)
        Q.append(Chromosome)
    return Q


result1 = graph_to_genome4([(2,4),(3,6),(5,1),(7,9),(10,12),(11,8)])
result2 = graph_to_genome4([(7,9),(10,12),(2,4),(3,6),(5,1),(11,8)])
result3 = graph_to_genome4([(5,1),(2,4),(3,6),(7,9),(10,12),(11,8)])
print(result1 == result2 == result3)  # False
print(graph_to_genome4([(2,4),(3,8),(7,5),(6,1)]))  # [['+3', '+1', '-2', '-4']]


def main():
    f = open('/home/wsl/rosalind/data/ba06i.txt', 'r')
    line = f.readline().strip('\n')
    temp = re.findall(r'\((.*?,.*?)\)', line)
    nodes = set()
    for frag in temp:
        fst, snd = [int(elem) for elem in frag.split(',')]
        nodes.add((fst,snd))
    f.close()

    start_time = time.time()
    answer = graph_to_genome(nodes)
    print("=====================================================")
    for ls in answer:
        print('(' + ' '.join(ls) + ')', end='')
    print()
    print(graph_to_genome({(8, 1), (4, 6), (5, 7), (3, 2)}))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
