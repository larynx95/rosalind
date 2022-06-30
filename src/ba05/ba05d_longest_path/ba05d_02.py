"""
Rosalind: BA5D
Find the Longest Path in a DAG

Longest Path in a DAG Problem
Find a longest path between two nodes in an edge-weighted DAG.

Given:
An integer representing the "source node" of a graph,
followed by an integer representing the "sink node" of the graph,
followed by an edge-weighted graph.
The graph is represented by a modified adjacency list in which the notation "0->1:7"
indicates that an edge connects node 0 to node 1 with weight 7.

Return:
The length of a longest path in the graph, followed by a longest path.
(If multiple longest paths exist, you may return any one.)

Sample Dataset
0
4
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3

Sample Output
9
0->2->3->4

═════════════════════════════════════════════════
References:
- solution by paps machine
"""

import time
import sys


def topsort(graph, indegree):
    def dfs(graph, start, visited, queue):
        visited[start] = True
        for v, w in graph[start]:
            if not visited[v]:
                dfs(graph, v, visited, queue)
        queue.insert(0,start)

    queue = []
    nodes = sorted(graph.keys())
    visited = {v:False for v in graph}
    while True:
        for v in nodes:
            if indegree[v] == 0 and not visited[v]:
                dfs(graph, start, visited, queue)
        else:
            break
    return queue


def LongestPath(DAG, queue, start, end):
    longest_path = {u: None for u in DAG.keys()}
    predecessors = {u: None for u in DAG.keys()}
    while len(queue) > 0:
        u = queue.pop(0)
        if longest_path[u] == None:
            longest_path[u] = 0
        for v, w in DAG[u]:
            if longest_path[v] == None or longest_path[u] + w > longest_path[v]:
                longest_path[v] = longest_path[u] + w
                predecessors[v] = u

    # construct the path using the predecessors table
    path = []
    u = end
    while u != -1:
        path.insert(0, u)
        if u == start:
            break
        u = predecessors[u]
    length = longest_path[end] - longest_path[start]
    return length, path


start = int(input())
end = int(input())
edges = [map(int, line.replace("->", ':').split(':')) for line in sys.stdin]

graph = {}
indegree = {}
for u, v, w in edges:
    if v not in indegree:
        indegree[v] = 1
    else:
        indegree[v] += 1
    if u not in indegree:
        indegree[u] = 0

    if u not in graph:
        graph[u] = [(v,w)]
    else:
        graph[u].append((v,w))
    if v not in graph:
        graph[v] = []

queue = topsort(graph, indegree)
length, path = LongestPath(graph, queue, start, end)
print(length)
print('->'.join(str(v) for v in path))