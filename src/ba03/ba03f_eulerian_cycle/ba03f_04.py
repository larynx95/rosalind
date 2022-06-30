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

https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA3F.py
'''

import sys
from re import split
from random import choice


def parse_adj_list(adj_list_text):
    adj_list = {}
    for elem in adj_list_text:
        temp = split(' -> ', elem)
        adj_list[temp[0]] = temp[1].split(',')
    return adj_list


def remove_edge(adj_list, from_node, to_node):
    adj_list[from_node].remove(to_node)
    if not adj_list[from_node]:
        del adj_list[from_node]
    return adj_list


def Eulerian_cycle(adj_list):
    # form a cycle Cycle by randomly walking in Graph
    start_node, edges = choice(list(adj_list.items()))
    target_node = choice(edges)
    adj_list = remove_edge(adj_list, start_node, target_node)

    Cycle = [start_node, target_node]
    current_node = target_node
    while current_node != start_node:
        edges = adj_list[current_node]
        target_node = choice(edges)
        adj_list = remove_edge(adj_list, current_node, target_node)
        current_node = target_node
        Cycle.append(current_node)

    while adj_list:
        potential_starts = [(idx, node) for idx, node in enumerate(Cycle) if node in adj_list]
        idx, new_start = choice(potential_starts)

        # form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking
        new_cycle = Cycle[idx:] + Cycle[1:idx + 1]

        target_node = choice(adj_list[new_start])
        adj_list = remove_edge(adj_list, new_start, target_node)
        current_node = target_node
        new_cycle.append(current_node)
        while current_node != new_start:
            edges = adj_list[current_node]
            target_node = choice(edges)
            adj_list = remove_edge(adj_list, current_node, target_node)
            current_node = target_node
            new_cycle.append(current_node)
        Cycle = new_cycle
    return Cycle


if __name__ == "__main__":
    '''
    Given: An Eulerian directed graph, in the form of an adjacency list.
    Return: An Eulerian cycle in this graph.
    '''
    input_lines = sys.stdin.read().splitlines()
    Adj_list = parse_adj_list(input_lines)

    print("->".join(Eulerian_cycle(Adj_list)))