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
- solution by Carlos Garcia
"""

import time


def LongestPath(Graph,source,sink):
    s = {}
    path = {}
    for i in range(len(Graph)):
        s[Graph[i][0]] = -float('inf')
        s[Graph[i][1]] = -float('inf')
    s[source] = 0
    for i in range(len(Graph)):
        if s[Graph[i][1]] < s[Graph[i][0]] + Graph[i][2]:
            s[Graph[i][1]] = s[Graph[i][0]] + Graph[i][2]
            path[Graph[i][1]] = Graph[i][0]
    Stack = [sink]
    Stack.append(path[sink])
    while Stack[-1] != source:
        Stack.append(path[Stack[-1]])
    return [s[sink],Stack[::-1]]


def main():

    f = open('/home/wsl/rosalind/data/ba05d.txt','r')
    a = f.read()
    a = a.split('\n')
    SOURCE = int(a[0])
    SINK = int(a[1])
    FROM = [int(n.split('->')[0]) for n in a[2:] if n != '']
    AUX = [n.split('->')[1] for n in a[2:] if n != '']
    TO = [int(n.split(':')[0]) for n in AUX]
    MARK = [int(n.split(':')[1]) for n in AUX]
    GRAPH = sorted(zip(FROM,TO,MARK))
    Solution = LongestPath(GRAPH,SOURCE,SINK)
    f.close()


    for i in range(len(Solution)):
        if i == 0:
            print(Solution[i])
        else:
            print('->'.join(map(str,Solution[i])))


if __name__ == "__main__":
    main()
