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

https://github.com/zhixue/rosalind_codes/blob/master/BA3F.py
'''

def parse_data(data, mmap):
	for line in data:
		temp = line.split(" -> ")
		mmap[temp[0]] = temp[1].split(",")
	return mmap

def combine_cycles(cycle1, cycle2):

	temp = cycle2.split("->")
	temp.pop()
	temp2 = cycle1.split("->")
	location = temp2.index(temp[0])
	temp2[location] = "*****"
	blah = "->".join(temp2)
	location = blah.find("*****")

	return cycle1[:location] + "->".join(temp) + "->" + cycle1[location:]

def find_cycle(cycles, mmap, key):
	cycle = ""
	if key == "none":
		# find node that has open edges and set it as start
		for key in mmap.keys():
			if mmap[key]:
				cycle = key
				break

		# set key to cycle

		key = cycle
	else:
		cycle = key

	# walk through until there is no edge left
	while mmap[key]:
		cycle += "->"
		# always choose first from the list
		key = mmap[key].pop(0)
		cycle += key
	if cycles != "":
		cycles = combine_cycles(cycles, cycle)
	else:
		cycles = cycle
	return (cycles, mmap)

def find_y(cycles, mmap):
	temp = cycles.split("->")
	for node in temp:
		if mmap[node] != []:
			return node
	return "completed"


with open('/home/wsl/rosalind/data/ba03f.txt') as f:
    ls = [x.rstrip() for x in f.readlines()]
    mmap = {}
    for one in ls:
        start = one.split(' ')[0]
        end = one.split(' ')[-1]
        if one.find(','):
            mmap[start] = end.split(',')
        else:
            mmap[start] = [end]



# cycles will hold all cycles found
cycles = ""


# first run through to find main cycle
temp = find_cycle(cycles, mmap, "none")
cycles = temp[0]
mmap = temp[1]
nextone = find_y(cycles, mmap)


while nextone != "completed":
    temp = find_cycle(cycles, mmap, nextone)
    cycles = temp[0]
    mmap = temp[1]
    nextone = find_y(cycles, mmap)
print(cycles)
