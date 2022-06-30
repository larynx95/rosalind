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

https://github.com/nutstick/rosalind/blob/master/ba3g.py
'''

g = {}
vertices = set()
indegree = {}
outdegree = {}

def join(l, sep):
    out = ''
    for element in l:
        out += str(element) + sep
    return out[:-len(sep)]

def euler_path(graph, vertices):
    # sum in-degree and out-degree
    indegree = dict.fromkeys(vertices, 0)
    outdegree = dict.fromkeys(vertices, 0)
    for vertex in vertices:
        if vertex in g:
            outdegree[vertex] = len(g[vertex])
            for adj in g[vertex]:
                indegree[adj] += 1

    # determine start and end point
    start = -1
    end = -1
    for vertex in vertices:
        if indegree[vertex] < outdegree[vertex]:
            start = vertex
        elif indegree[vertex] > outdegree[vertex]:
            end = vertex
    print(start, end)

    # tour
    current_path, circuit, v = [start], [], start
    while len(current_path) > 0:
        if outdegree[v]:
            current_path.append(v)
            nextv = g[v].pop()
            outdegree[v] -= 1
            v = nextv
        else:
            circuit.append(v)
            v = current_path.pop()
    # print(current_path, circuit)
    circuit.reverse()
    return circuit

with open('/home/wsl/rosalind/data/ba03g_output.txt', 'w') as out:
    with open('/home/wsl/rosalind/data/ba03g.txt', 'r') as f:
        for line in f:
            nodes = line[:-1].split(' -> ')
            u = int(nodes[0])
            vs = [int(n) for n in nodes[1].split(',')]
            vertices.add(u)
            vertices = vertices  | set(vs)
            g[u] = vs
        path = euler_path(g, list(vertices))
        out.write(join(path, '->'))