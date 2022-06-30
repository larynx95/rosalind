"""
Rosalind: BA7B
Compute Limb Lengths in a Tree

Limb Length Problem
Find the limb length for a leaf in a tree.

Given:
An integer n, followed by an integer j between 0 and n-1,
followed by a space-separated additive distance matrix D
(whose elements are integers).

Return:
The limb length of the leaf in Tree(D)
corresponding to row j of this distance matrix (use 0-based indexing).

Sample Dataset
4
1
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0

Sample Output
2

═════════════════════════════════════════════════

Info.
  * Tree(D) --> distance matrix (D) --> limb length
      1      2                0   1   2   3
     2↕  4   ↕6   BA7A     -----------------   BA7B    i \     / ?
      4 <--> 5    ====>    0  0   13  21  22   ===>       ?...?
    11↕      ↕7            1  13  0   12  13           ? /     \ ?
      0      3             2  21  12  0   13           TODO: How to get limb length without any info about internal node?
                           3  22  13  13  0
  - TODO: Questions:
    - Does neighbors always have minimum length?
      If D(i,j) has minimum value in distance matrix D,
      are leaves i and j always neighbors? No. Look at the example above.
      The min value is 12 in distance matrix (D(1,2) or D(2,1)), but 1 and 2 are not neighbors.
    - Can I get only one limb length from distance matrix D?
      Or do I need to construct whole simple tree before getting a specific limb length?
  - limb length = distance from a leaf i to Parent(i)
  - what I can get from given data
    - There're max and min values in the given distance matrix (D).
      The minimum D(i,j) doesn't always mean that leaves i and j are neighbors.
      The maximum D(i,j) always means that leaves i and j can't be neighbors to each other.
    - in D(i,j), i != j
    - Theorem: Every simple tree at least three nodes has a pair of neighboring leaves.

  ┌─────────────────────────────────────────────────────────────┐
  │ Distance-Based Phylogeny Problem:                           │
  │ Reconstruct an evolutionary tree fitting a distance matrix. │
  │   Input: A distance matrix.                                 │
  │   Output: A tree fitting this distance matrix.              │
  └─────────────────────────────────────────────────────────────┘

  * how to get limb length
    - example (limb length of leave i)
      A. select minimum length of D(i,j) = i, j are neighbers
         TODO: Does neighbors always have minimum length?
               Are leaves i and j always neighbors?

       i \       / k    Any leave can be k.
          m --- ?       D(i,j) is minimum value in distance matrix D.
       j /       \ ?    conditions:
                          - leaves i, j sharing a parent node (neighboring leaves)

     d(k,m) = [(d(i,m) + d(k,m)) + (d(j,m) + d(k,m)) - (d(i,m) + d(j,m))] / 2
     d(k,m) = [d(i,k) + d(j,k) - d(i,j)] / 2
     d(k,m) = [D(i,k) + D(j,k) - D(i,j)] / 2
     d(i,m) = D(i,k) - [D(i,k) + D(j,k) - D(i,j)] / 2
     d(i,m) = [D(i,k) + D(i,j) - D(j,k)] / 2  <-- This is limb length of leave i.

     B. select not minimum length of D(i,j) == i, j are not neighbors
        TODO: What if i and j are not neighbors?

       i \       / j    Any leave can be k.
          m --- ?       D(i,j) is NOT minimum length in distance matrix D.
       k /       \ ?

     d(j,m) = [(d(i,m) + d(j,m)) + (d(k,m) + d(j,m)) - (d(i,m) + d(k,m))] / 2
     d(j,m) = [d(i,j) + d(j,k) - d(i,k)] / 2
     d(j,m) = [D(i,j) + D(j,k) - D(i,k)] / 2
     d(i,m) = D(i,j) - D(j,m) = D(i,j) - [D(i,j) + D(j,k) - D(i,k)] / 2
     d(i,m) = [D(i,j) + D(i,k) - D(j,k)] / 2  <-- This is limb length of leave i.
     d(j,?) = ??  <-- I can't get the limb length of j in this way.

  * TODO: construct simple Tree (T) from distance matrix (D)
    - distance matrix (D) --> simple Tree (T)
    - toward a recursive algorithm (textbook p.12)
    - steps:
       a. find a pair of neighboring leaves i and j
          by selecting the "minimum" element D(i,j) in the distance matrix;
          TODO: Is this always right?
       b. replace i and j with their parent,
          and recompute the distances from this parent to all other leaves as described above;
       c. solve the Distance-Based Phylogeny problem for the smaller tree;
       d. add the previously removed leaves i and j back to the tree.

Plan 1.
  * using equaltion in textbook (p.12~13)
    for neighboring leaves i and j sharing a parent node m,
    the following equality holds for every other leaf k in the tree:

    d(k,m) = [(d(i,m)+d(k,m)) + (d(j,m)+d(k,m)) - (d(i,m)+d(j,m))] / 2
           = [d(i,k) + d(j,k) - d(i,j)] / 2

    But i and j are not always neighbors. Be careful.

  * find neighboring two leaves --> TODO: Is this easy?

Plan 2.
  * focusing on MIN value
  * It is not easy to find neighboring two leaves from the given distance matrix D.
    Leaves i and j can be neighbors or not.
    But if leaves i and j are neighbors, the value of d(i,j) will be minimum.

  * get all limb length and select min value
    For leaf i, limb length of leave i must be the smallest value among d(i,_).

Plan 3.
  * focusing on MAX value
  * Even if D(i,j) is minimum value in distance matrix D,
    it doesn't always mean that leaves i and j are neighbors.
    But if D(i,j) is maximum value, it means that leaves i and j are not neighbors.

═════════════════════════════════════════════════

References:
- Toward an Algorithm for Distance-Based Phylogeny Construction
  https://www.youtube.com/watch?v=TSIkxBoaykA
- Split a string with unknown number of spaces as separator in Python
  https://stackoverflow.com/questions/4309684/split-a-string-with-unknown-number-of-spaces-as-separator-in-python
"""
#!/usr/bin/env python
import time
import random
import re


def limb_length_wrong(matrix, leaf_i):
    """
    ([[int]],int) -> int
    returns a limb length for a leaf
    """
    # get two leaves from a distance matrix
    n = len(matrix)
    leaf_j = random.choice([num for num in range(n) if num != leaf_i])
    leaf_k = random.choice([num for num in range(n) if num != leaf_i and num != leaf_j])
    # get three distances - d(i,j), d(i,k), d(k,j)
    d_ij = matrix[leaf_i][leaf_j]
    d_ik = matrix[leaf_i][leaf_k]
    d_jk = matrix[leaf_j][leaf_k]
    # get limb length [D(i,k) + D(i,j) - D(j,k)] / 2
    answer = (d_ik + d_ij - d_jk) // 2
    return answer


def limb_length(matrix, i):
    """
    ([[int]],int) -> int
    returns a limb length for a leaf
    focusing on minimum value in distance matrix D
    """
    n = len(matrix)
    min_limb_length = float('inf')
    for j in range(n):
        for k in range(n):
            if i != j and j != k and i != k:
                limb_length = (matrix[i][k] + matrix[i][j] - matrix[j][k]) // 2
                if limb_length < min_limb_length:
                    min_limb_length = limb_length
    return min_limb_length


def limb_length2(matrix, j):
    """
    ([[int]],int) -> int
    returns a limb length for a leaf
    list comprehension version
    """
    n = len(matrix)
    return min(matrix[i][j] + matrix[k][j] - matrix[i][k] for i in range(n) for k in range(n) if i != j and k != j) // 2


def limb_length3(matrix, leaf):
    """
    ([[int]],int) -> int
    returns a limb length for a leaf
    focusing on maximum value in distance matrix D
    """
    n = len(matrix)
    min_limb_length = max(matrix[0])  # starting from max value in distance matrix D
    for j in range(n):
        # diagonal elements
        if j == leaf:
            continue
        for k in range(n):
            if k == leaf or j == k:
                continue
            c = (matrix[leaf][j] + matrix[leaf][k] - matrix[j][k]) / 2
            if c < min_limb_length:
                min_limb_length = c
    return int(min_limb_length)


def main():
    f = open('/home/wsl/rosalind/data/ba07b.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    n = int(lines[0])
    leaf_i = int(lines[1])
    matrix = []
    for line in lines[2:]:
        row = [int(elem) for elem in line.split()]
        matrix.append(row)
    f.close()

    start_time = time.time()
    print(limb_length(matrix, leaf_i))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()