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
- solution by other person
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


def genome_to_colored(genome):
    """
    [[int]] -> {(int,int)}
    returns a list of list of integers
    >>> genome_to_colored([[1,-2,-3],[4,5,-6]])
        {(2,4),(3,6),(5,1),(8,9),(10,12),(11,7)}
    """
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(len(chromosome)):
            edges.add((nodes[2*j+1], nodes[2*(j+1) % len(nodes)]))
    return edges


def two_break_distance(p, q):
    """
    solution by other person
    """
    edges = genome_to_colored(p).union(genome_to_colored(q))
    blocks = sum(map(len, p))
    graph = {_+1: [] for _ in range(2 * blocks)}
    for edge in edges:
        u, v = edge
        graph[u].append(v)
        graph[v].append(u)
    unvisited = [1] * 2 * blocks
    # number of cycles is the same as the number of connected components for P + Q graph
    cnt = 0
    while any(unvisited):
        queue = [unvisited.index(1) + 1]
        unvisited[queue[0] - 1] = 0
        cnt += 1
        while queue:
            u = queue.pop()
            for v in graph[u]:
                if unvisited[v-1]:
                    unvisited[v-1] = 0
                    queue.append(v)
    return blocks - cnt


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
    f = open('/home/wsl/rosalind/data/ba06c.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    genome_p = read_genome(lines[0])
    genome_q = read_genome(lines[1])
    f.close()

    start_time = time.time()
    # print(two_break_distance([[+1,+2,+3,+4,+5,+6]],[[+1,-3,-6,-5],[+2,-4]]))
    print(two_break_distance(genome_p, genome_q))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()