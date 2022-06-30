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

https://github.com/raomanus/rosalind/blob/master/ba3g.py
'''

def constructGraph(edges):
	graph = dict()
	inDegree = dict()
	outDegree = dict()
	for edge in edges:
		s,e = edge.split('->')
		s = s.strip()
		e = e.strip().split(',')
		try:
			graph[s].extend(e)
		except:
			graph[s] = e
		try:
			outDegree[s] += len(e)
		except:
			outDegree[s] = len(e)
		for end in e:
			try:
				inDegree[end] += 1
			except:
				inDegree[end] = 1
	return(graph,outDegree,inDegree)

def findEulerPath(graph, start):
	path = []
	stack = [start]
	location = graph[start].pop()
	keys = list(graph.keys())
	while len(stack) != 0:
		if (location not in keys) or (len(graph[location]) == 0):
			path.append(location)
			location = stack.pop()
		else:
			stack.append(location)
			location = graph[location].pop()
	path.append(start)
	return path[::-1]

file = open("/home/wsl/rosalind/data/ba03g.txt","r")
edges = file.read().splitlines()
start = None

graph,outDegree,inDegree = constructGraph(edges)
inKeys = list(inDegree.keys())
outKeys = list(outDegree.keys())

for ik in outKeys:
	if ik in inKeys:
		if inDegree[ik] < outDegree[ik]:
			start = ik
	else:
		start = ik

path = findEulerPath(graph, start)
output = path[0]

for node in path[1:]:
	output += "->"+node

print(output)