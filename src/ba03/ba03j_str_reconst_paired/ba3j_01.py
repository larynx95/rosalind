"""
Rosalind: BA3J
Reconstruct a String from its Paired Composition

Since increasing read length presents a difficult experimental problem,
biologists have suggested an indirect way of increasing read length by generating read-pairs,
which are pairs of reads separated by a fixed distance d in the genome.

You can think about a read-pair as a long "gapped" read of length k + d + k
whose first and last k-mers are known but whose middle segment of length d is unknown.
Nevertheless, read-pairs contain more information than k-mers alone,
and so we should be able to use them to improve our assemblies.
If only you could infer the nucleotides in the middle segment of a read-pair,
you would immediately increase the read length from k to 2 · k + d.

Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d.
We use the notation (Pattern1|Pattern2) to refer to a a (k,d)-mer whose k-mers are Pattern1 and Pattern2.
The (k,d)-mer composition of Text, denoted PairedCompositionk,d(Text),
is the collection of all (k,d)- mers in Text (including repeated (k,d)-mers).

String Reconstruction from Read-Pairs Problem
Reconstruct a string from its paired composition.

Given: Integers k and d followed by a collection of paired k-mers PairedReads.

Return: A string Text with (k, d)-mer composition equal to PairedReads.
        (If multiple answers exist, you may return any one.)

Sample Dataset
4 2
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT

Sample Output
GTGGTCGTGAGATGTTGA

═════════════════════════════════════════════════

"""

import time
import sys


def find_start(al):
    for u in al:
        out_degree = len(al[u])
        in_degree = 0
        for v in al:
            in_degree += al[v].count(u)
        if in_degree + out_degree % 2 == 1:
            return u


def are_unexplored_edges(ue):
    flag = False
    for u in ue:
        if ue[u]:
            flag = True
    return flag


def random_path_walk(start, al, ue):
    # path walk
    path = [start]
    while True:
        u = path[-1]
        if u not in al:
            break
        for i in range(len(al[u])):
            if i in ue[u]:
                v = al[u][i]
                if v == start:
                    continue
                ue[u].remove(i)
                path.append(v)
                break
        if path[-1] == u:
            break
    return path


def random_cycle_walk(start, al, ue):
    # cyclic walk
    cycle = [start]
    while True:
        u = cycle[-1]
        if u not in al:
            break
        for i in range(len(al[u])):
            if i in ue[u]:
                v = al[u][i]
                ue[u].remove(i)
                cycle.append(v)
                break
        if cycle[-1] == start:
            break
    return cycle


def eulerian_path(adj_list):
    start = find_start(adj_list)
    unexplored_edges = {u: set(range(len(adj_list[u]))) for u in adj_list}
    # initial path
    path = random_path_walk(start, adj_list, unexplored_edges)
    while are_unexplored_edges(unexplored_edges):
        for i in range(len(path)):
            u = path[i]
            if u not in adj_list:
                continue
            if unexplored_edges[u]:
                cycle = random_cycle_walk(u, adj_list, unexplored_edges)
                path = path[:i] + cycle + path[i+1:]
    return path


def paired_string_reconstruction(pc, k, d):
    graph = {}
    for pair in pc:
        prefix = (pair[0][:k-1], pair[1][:k-1])
        suffix = (pair[0][1:], pair[1][1:])
        if prefix in graph:
            graph[prefix].append(suffix)
        else:
            graph[prefix] = [suffix]
    path = eulerian_path(graph)
    string1, string2 = '', ''
    for pair in path:
        read1, read2 = pair
        string1 += read1[-1] if string1 else read1[:(k-1)]
        string2 += read2[-1] if string2 else read2[:(k-1)]
    string = string1 + string2[-(k+d):]
    return string


def main():
    f = open('/home/wsl/rosalind/data/ba03j.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, d = [int(num.strip()) for num in lines[0].split()]
    paired_reads = []
    for line in lines[1:]:
        prefix, suffix = [rd.strip() for rd in line.split('|')]
        paired_reads.append((prefix, suffix))
    f.close()

    start_time = time.time()
    print(paired_string_reconstruction(paired_reads,k,d))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()