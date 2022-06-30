'''
Rosalind: BA3F
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

    EULERIANCYCLE(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle' by traversing Cycle (starting at newStart) and then randomly walking
            Cycle <- Cycle'
        return Cycle

Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An Eulerian directed graph, in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3
1 -> 0
2 -> 1,6
3 -> 2
4 -> 2
5 -> 4
6 -> 5,8
7 -> 9
8 -> 7
9 -> 6

Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

-------------------------------------------------

https://github.com/Leoberium/BA/blob/master/Chapter3/BA3F.py
'''

import sys


def are_unexplored_edges(ue):
    flag = False
    for u in ue:
        if ue[u]:
            flag = True
    return flag


def random_walk(start, al, ue):
    # cyclic walk
    cycle = [start]
    v = start
    while True:
        u = cycle[-1]
        for i in range(len(al[u])):
            if i in ue[u]:
                v = al[u][i]
                ue[u].remove(i)
                cycle.append(v)
                break
        if v == start:
            break
    return cycle


def eulerian_cycle(adj_list):
    start = min(adj_list.keys())
    unexplored_edges = {u: set(range(len(adj_list[u]))) for u in adj_list}
    # initial cycle
    cycle = random_walk(start, adj_list, unexplored_edges)
    while are_unexplored_edges(unexplored_edges):
        for i in range(len(cycle)):
            u = cycle[i]
            if unexplored_edges[u]:
                start = u
                new_cycle = random_walk(start, adj_list, unexplored_edges)
                cycle = cycle[:i] + new_cycle + cycle[i+1:]
    return cycle


def main():
    f = open('/home/wsl/rosalind/data/ba03f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    adj_list = dict()
    # for line in sys.stdin:
    for line in lines:
        line = line.strip()
        tail, head = line.split(' -> ')
        adj_list[int(tail)] = list(map(int, head.split(',')))
    cycle = eulerian_cycle(adj_list)
    print('->'.join(map(str, cycle)))


if __name__ == '__main__':
    main()