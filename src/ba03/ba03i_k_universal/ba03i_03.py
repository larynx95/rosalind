"""
Rosalind: BA3I
Find a k-Universal Circular String

A k-universal circular string is a circular string
that contains every possible k-mer constructed over a given alphabet.

k-Universal Circular String Problem
Find a k-universal circular binary string.

Given: An integer k.

Return: A k-universal circular string.
(If multiple answers exist, you may return any one.)

Sample Dataset
4

Sample Output
0000110010111101

═════════════════════════════════════════════════

https://www.geeksforgeeks.org/de-bruijn-sequence-set-1/
"""
#!/usr/bin/env python
import time


# Python3 implementation of
# the above approach
import math

seen = set()
edges = []

# Modified DFS in which no edge is traversed twice
def dfs( node, k, A):

    for i in range(k):
        str = node + A[i]
        if (str not in seen):
            seen.add(str)
            dfs(str[1:], k, A)
            edges.append(i)

# Function to find a de Bruijn sequence of order n on k characters
def deBruijn(n, k, A):

    # Clearing global variables
    seen.clear()
    edges.clear()

    startingNode = A[0] * (n - 1)
    dfs(startingNode, k, A)

    S = ""

    # Number of edges
    l = int(math.pow(k, n))
    for i in range(l):
        S += A[edges[i]]

    S += startingNode
    return S

# Driver code
n = 3
k = 2
A = "01"

print(deBruijn(n, k, A))

# This code is contributed by shubhamsingh10

def main():
    f = open('/home/wsl/rosalind/data/ba03i.txt', 'r')
    k = int(f.readline().strip())
    f.close()

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
