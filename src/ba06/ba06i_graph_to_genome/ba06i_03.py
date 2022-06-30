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
def read_data(file_name):
    with open(file_name, "r") as file:
        line = file.readline().strip().split(", (")

    length = len(line)
    line[0] = line[0].replace("(", "")
    graph = []
    for i in range(length):
        line[i] = line[i].replace(")", "")
        node1, node2 = line[i].split(",")
        graph.append(int(node1))
        graph.append(int(node2))

    return graph


def cycle_to_chrom(nodes):
    length = len(nodes)
    chrom = [0]*(length // 2)
    for j in range(length // 2):
        if nodes[2*j] < nodes[2*j+1]:
            chrom[j] = nodes[2*j+1]/2
        else:
            chrom[j] = -nodes[2*j]/2

    return chrom


def _get_pair_elem(start):
    if start % 2 == 0:
        return start-1
    else:
        return start+1


def graph_to_genome(graph):
    genom = []
    start = graph[0]
    st_index = 0
    end = _get_pair_elem(start)
    for i in range(len(graph) // 2):
        if graph[2*i+1] == end:
            cycle = [end] + graph[st_index:2*i+2]
            genom.append(cycle_to_chrom(cycle))
            if 2*i+2 < len(graph):
                start = graph[2*i+2]
                st_index = 2*i+2
                end = _get_pair_elem(start)
    return genom


def form_answ(genom):
    answ = ""
    for j in range(len(genom)):
        for i in range(len(genom[j])):
            if genom[j][i] > 0:
                genom[j][i] = "+" + str(genom[j][i])
            else:
                genom[j][i] = str(genom[j][i])

        answ += "(" + " ".join(genom[j]) + ")"

    print (answ)


def main():
    data = read_data("/home/wsl/rosalind/data/ba06i.txt")
    genom = graph_to_genome(data)
    form_answ(genom)


if __name__ == '__main__':
    main()
