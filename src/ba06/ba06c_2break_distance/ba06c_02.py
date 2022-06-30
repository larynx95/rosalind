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
* solution by other person 'root'

* Specifically, the 2-Break Distance Theorem says d(P,Q)=blocks(P,Q)−cycles(P,Q).
  Getting blocks(P,Q) is easy enough, as it is just the total number of elements in P (or Q).
  For cycles(P,Q), note that each connected component of the breakpoint graph is a unique cycle,
  so the number of cycles is just the number of connected components,
  and hence cycles(P,Q)=#Connect Components.
"""
#!/usr/bin/env python
import time
from collections import defaultdict

def two_break_dist(P, Q):
    """
    ([[int]],[[int]]) -> int
    """
    # Construct the break point graph of P and Q.
    graph = defaultdict(list)
    for perm_cycle in P+Q:
        n = len(perm_cycle)
        for i in range(n):
            # Add the edge between consecutive items (both orders since the breakpoint graph is undirected).
            # Note: Modulo n in the higher index for the edge between the last and first elements.
            graph[perm_cycle[i]].append(-1*perm_cycle[(i+1) % n])
            graph[-1*perm_cycle[(i+1) % n]].append(perm_cycle[i])
    print(graph)
    # BFS to find the number of connected components in the breakpoint graph.
    component_count = 0
    remaining = set(graph.keys())
    while remaining:
        component_count += 1
        queue = {remaining.pop()}
        while queue:
            current = queue.pop()
            new_nodes = {node for node in graph[current] if node in remaining}
            queue |= new_nodes
            remaining -= new_nodes
    # Theorem: d(P,Q) = blocks(P,Q) - cycles(P,Q)
    return sum(map(len,P)) - component_count


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
    print(two_break_dist(genome_p, genome_q))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()