"""
get all paths from a graph

Example:
  graph = {6:[7,15], 15:[9,16], 16:[21], 9:[13], 13:[4], 4:[1], 1:[5]}


Output:
  6->7
  6->15->9->13->4->1->5
  6->15->16->21

═════════════════════════════════════════════════

References:
- Find all paths without specifing end node?
  https://stackoverflow.com/questions/33196046/find-all-paths-without-specifing-end-node
"""

import time


def allpaths1(graph, start, path=[]):
    path = path + [start]
    if start not in graph:
        return [path]
    paths = [path]
    for node in graph[start]:
        if node not in path:
            newpaths = allpaths1(graph, node, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


graph = {6:[7,15], 15:[9,16], 16:[21], 9:[13], 13:[4], 4:[1], 1:[5]}
print(allpaths1(graph, 6))


def allpaths2(graph, start, path=[]):
    path = path + [start]
    yield path
    if start not in graph:
        return
    for node in graph[start]:
        if node not in path:
            yield from allpaths2(graph, node, path)


graph = {6:[7,15], 15:[9,16], 16:[21], 9:[13], 13:[4], 4:[1], 1:[5]}
print(list(allpaths2(graph, 6)))


def main():
    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
