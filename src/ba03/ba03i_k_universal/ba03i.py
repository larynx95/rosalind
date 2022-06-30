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

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> Construct De Bruijn Graph (BA3D)
      -> Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> Eulerian Cycle (BA3F)
        -> HERE: k-universal string problem (BA3I)
      -> PREV: Eulerian Path (BA3G)
      ↓
    * Read-Pair (BA3J)
      -> NEXT: Construct a String Spelled by a Gapped Genome Path (BA3L)


Plan 1.
* steps:
  i.   create a list k-mers - cartesian product with "01"
       All (not some) of these binary k-mers should be used to reconstruct k-universal string.
  ii.  get De Bruijn graph from the list (exclude edges to node itself)
       split each k-mer into prefix and suffix (k-1)-mers
       then construct De Bruijn graph represented as dictionary
  iii. get Eulerian cycle from De Bruijn graph
  iv.  find out the way how to reconstruct string from Eulerian cycle
       BA3B string reconstruction function is not for circular. Be careful!
       ['00','01','11','11','10','01','10','00','00']

       index     0      1 ...     -3      -2      -1
       ---------------------------------------------
       k=2       0      1          1       0   /   0
       k=3      00     01 ...     10   /  00      00
       k=4     000    001 ...   10/0     000     000
       k=5    0000   0001 ...  10/00    0000    0000
       k=6   00000  00001 ... 10/000   00000   00000
                              -- ---
                             (a)  (b)   <-- (a) is always '10', but (b) changes

═════════════════════════════════════════════════

References:
- fast factorial problem
  - Why python math.factorial(x) is very fast?
    https://stackoverflow.com/questions/11313998/why-python-math-factorialx-is-very-fast
  - Tabulation vs Memoization
    https://www.geeksforgeeks.org/tabulation-vs-memoization/
  - Dynamic Programming Tutorial: making efficient programs in Python
    https://www.educative.io/blog/python-dynamic-programming-tutorial
- Generate all permutation of a set in Python
  https://www.geeksforgeeks.org/generate-all-the-permutation-of-a-list-in-python/
- How to generate all permutations of a list?
  https://stackoverflow.com/questions/104420/how-to-generate-all-permutations-of-a-list
- How to print like printf in Python3?
  https://stackoverflow.com/questions/19457227/how-to-print-like-printf-in-python3
"""
#!/usr/bin/env python
import time
import itertools


def prod_itertools(text, k):
    """
    (str,int) -> [str]
    >>> prod_itertools("ACG",2)
        ['AA','AC','AG','CA','CC','CG','GA','GC','GG']
    """
    from itertools import product
    result = []
    for tup in itertools.product(text, repeat=k):
        temp = ''
        for elm in tup:
            temp += elm
        result.append(temp)
    return result


def prod_recur(text, k):
    """
    (str,int) -> gen
    >>> list(prod_recur("ACG",2))
        ['AA','AC','AG','CA','CC','CG','GA','GC','GG']
    """
    if k == 1:
        for ch in text:
            yield ch
    else:
        for ch in text:
            for s in prod_recur(text, k - 1):
                yield ch + s


def prod_lscomp(text, k):
    """
    (str,int) -> [str]
    >>> prod_lscomp('01',3)
        ['000','001','010','011','100','101','110','111']
    """
    return ["".join(a) for a in itertools.product(text, repeat = k)]


def de_bruijn(kmers):
    """
    [a] -> {a:[a]}
    returns DeBruijnk graph from k-mer patterns (BA3E)
    representation of De Bruijn graph as list of tuples: not good idea, use dictionary
    >>> de_bruijn(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])
        {'GAG':['AGG'],'CAG':['AGG','AGG'],'GGG':['GGG','GGA'],'AGG':['GGG'],'GGA':['GAG']}
    >>> de_bruijn( ['000','001','010','011','100','101','110','111'])
        {'00':['00','01'],'01':['10','11'],'10':['00','01'],'11':['10','11']}
    """
    result = {}
    for kmer in kmers:
        key = kmer[:-1]
        val = kmer[1:]
        if key not in result:
            result[key] = [val]
        else:
            result[key].append(val)
    return result


def eulerian_cycle(graph):
    """
    {a:[a]} -> [a]
    returns Eulerian cycle
    """
    cycle = [min(graph.keys())]
    while len(graph) > 0:
        if cycle[0] == cycle[-1]:
            while not cycle[0] in graph:
                cycle.pop(0)
                cycle.append(cycle[0])
        source = cycle[-1]
        cycle.append(graph[source].pop())
        if len(graph[source]) == 0: del graph[source]
    return cycle


def str_reconst(kmers):
    """
    [str] -> str
    returns a merged string from overlapping k-mers (BA3B)
    >>> str_reconst(['00','01','11','11','10','01','10','00','00'])
        00011101
    """
    result = kmers[0]
    for node in kmers[1:]:
        result += node[-1]
    return result


def k_universal_strings(k):
    """
    int -> {str}
    returns k-universal circular string
    """
    kmers = prod_lscomp('01',k)
    graph = de_bruijn(kmers)
    cycle = eulerian_cycle(graph)
    indices = [i for i, x in enumerate(cycle) if x == cycle[0]]
    result = ""
    if k > 3:
        result = str_reconst(cycle[:-2])[:-(k-3)]
    return result


def main():
    f = open('/home/wsl/rosalind/data/ba03i.txt', 'r')
    k = int(f.readline().strip())
    f.close()

    start_time = time.time()
    print(k_universal_strings(k))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
