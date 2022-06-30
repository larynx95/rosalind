'''
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node is called an Eulerian path.

Eulerian Path Problem
Find an Eulerian path in a graph.

Given: A directed graph that contains an Eulerian path,
       where the graph is given in the form of an adjacency list.

Return: An Eulerian path in this graph.

Sample Dataset
0 -> 2
1 -> 3
2 -> 1
3 -> 0,4
6 -> 3,7
7 -> 8
8 -> 9
9 -> 6

Sample Output
6->7->8->9->6->3->0->2->1->3->4

-------------------------------------------------

https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA3G.py
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

def Eulerian_path(adj_list):
    deg_diffs = {}
    for source, targets in adj_list.items():
        if source in deg_diffs:
            deg_diffs[source] += len(targets)
        else:
            deg_diffs[source] = len(targets)
        for target in targets:
            if target in deg_diffs:
                deg_diffs[target] -= 1
            else:
                deg_diffs[target] = -1

    to_add_s = [node for node, diff in deg_diffs.items() if diff == -1][0]
    to_add_t = [node for node, diff in deg_diffs.items() if diff == 1][0]
    if to_add_s in adj_list:
        adj_list[to_add_s].append(to_add_t)
    else:
        adj_list[to_add_s] = [to_add_t]
    cycle = Eulerian_cycle(adj_list)
    idx = 0
    while True:
        if cycle[idx] == to_add_s and cycle[idx + 1] == to_add_t:
            break
        idx += 1
    return cycle[idx + 1:] + cycle[1:idx + 1]


if __name__ == "__main__":
    '''
    Given: A directed graph that contains an Eulerian path, where the graph is given in the form of an adjacency list.
    Return: An Eulerian path in this graph.
    '''
    # input_lines = sys.stdin.read().splitlines()
    f = open('/home/wsl/rosalind/data/ba03g.txt', 'r')
    input_lines = [line.strip() for line in f.readlines()]
    gdict = parse_adj_list(input_lines)
    f.close()
    print("->".join(Eulerian_path(gdict)))