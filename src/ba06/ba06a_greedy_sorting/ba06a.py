"""
Rosalind: BA6A
Implement GreedySorting to Sort a Permutation by Reversals

Implement GreedySorting
Given:
A signed permutation P.

Return:
The sequence of permutations corresponding to applying GreedySorting to P,
ending with the "identity permutation".

Sample Dataset
(-3 +4 +1 +5 -2)

Sample Output
(-1 -4 +3 +5 -2)
(+1 -4 +3 +5 -2)
(+1 +2 -5 -3 +4)
(+1 +2 +3 +5 +4)
(+1 +2 +3 -4 -5)
(+1 +2 +3 +4 -5)
(+1 +2 +3 +4 +5)

═════════════════════════════════════════════════

    [ Where am I? ]

    * HERE: GreedySort genes in a chromosome (BA6A)
      ↓
    * Helper function for 2-Break sorting algorithm
      -> NEXT: ChromosomeToCycle (BA6F)

Info.
  * identity permutation:
    ; ordered from smallest to largest with positive directions

  * terminology
                                         genome :: [[signed int]]
                                     chromosome :: [signed int]
    synteny blocks, black edges, directed edges :: {(int,int)}
                colored edges, undirected edges :: {(int,int)}
    - k-sorted: first k-1 elements are sored, but element k is unsorted
    - k-sorting reversal: for every (k-1)-sorted permutation P,
      there exists a single reversal, called the k-sorting reversal,
      that fixes the first k-1 elements of P and moves element k to the k-th position
    - d_rev(P): the minimum number of reversals required to sort P into the identity permutation as d_rev(P)

  * sample dataset
    [-3 +4 +1]+5 -2
    [-1]-4 +3 +5 -2
     +1[-4 +3 +5 -2]
     +1 +2[-5 -3]+4
     +1 +2 +3[+5 +4]
     +1 +2 +3[-4]-5
     +1 +2 +3 +4[-5]
     +1 +2 +3 +4 +5

  * another example
    +1 [-7  +6 -10  +9  -8  +2]-11  -3  +5  +4      reversal
    +1 [-2] +8  -9 +10  -6  +7 -11  -3  +5  +4      (-) to (+)
    +1  +2 [+8  -9 +10  -6  +7 -11  -3] +5  +4      reversal
    +1  +2  +3[+11  -7  +6 -10  +9  -8  +5  +4]     reversal
    +1  +2  +3 [-4] -5  +8  -9 +10  -6  +7 -11      (-) to (+)
    +1  +2  +3  +4 [-5] +8  -9 +10  -6  +7 -11      (-) to (+)
    +1  +2  +3  +4  +5 [+8  -9 +10  -6] +7 -11      reversal
    +1  +2  +3  +4  +5  +6[-10  +9  -8  +7]-11      reversal
    +1  +2  +3  +4  +5  +6 [-7] +8  -9 +10 -11      (-) to (+)
    +1  +2  +3  +4  +5  +6  +7  +8 [-9]+10 -11      (-) to (+)
    +1  +2  +3  +4  +5  +6  +7  +8  +9 +10[-11]     (-) to (+)
    +1  +2  +3  +4  +5  +6  +7  +8  +9 +10 +11

Plan 1.
  * sort abs(val) in list
    -3 +4 +1 +5 -2   => 1  2  3  4  5    (all elements are unique)
    TODO: Study other sorting algorithms!

  * compare, set interval for reversal --> erversal, sign change
      0  1  2  3  4   <-- zero based index
      1  2  3  4  5   idx  match  idx_r  interval  action
    ----------------------------------------------------
    [-3]+4 +1 +5 -2   0    N      2      [0:2+1]   reversal
    [-1]-4 +3 +5 -2   0    Y      -      -         change direction (sign)
     +1[-4 +3 +5 -2]  1    N      4      [1:4+1]   reversal
     +1 +2[-5 -3]+4   2    N      3      [2:3+1]   reversal
     +1 +2 +3[+5 +4]  3    N      4      [3:4+1]   reversal
     +1 +2 +3[-4]-5   3    Y      -      -         change direction (sign)
     +1 +2 +3 +4[-5]  4    Y      -      -         change direction (sign)
     +1 +2 +3 +4 +5   completed!

    idx  incorrect => reverse, change sign
    sign incorrect => change sign

  * helper functions
    - sorting positive numbers: required and consider efficient algorithm
      TODO: What is the most efficient algorithm for sorting?
    - string reversal, changing sign: not required - use `[::-1]` syntax

Plan 2.
  * algorithm
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║ GREEDYSORTING(P)                                                      ║
  ║     approxReversalDistance <- 0                                       ║
  ║     for k <- 1 to |P|                                                 ║
  ║         if element k is not sorted                                    ║
  ║             apply the k-sorting reversal to P                         ║
  ║             approxReversalDistance <- approxReversalDistance + 1      ║
  ║             if the k-th element of P is -k                            ║
  ║                 apply the k-sorting reversal to P                     ║
  ║                 approxReversalDistance <- approxReversalDistance + 1  ║
  ║     return approxReversalDistance                                     ║
  ╚═══════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- chapter 05. merge, mergesort
╔═════════════════════════════════════════════════════════════════════════════════════╗
║ MERGE(List1, List2)                                                                 ║
║   SortedList empty list                                                             ║
║   while both List1 and List2 are non-empty                                          ║
║     if the smallest element in List1 is smaller than the smallest element in List2  ║
║       move the smallest element from List1 to the end of SortedList                 ║
║     else                                                                            ║
║       move the smallest element from List2 to the end of SortedList                 ║
║   move any remaining elements from either List1 or List2 to the end of SortedList   ║
║   return SortedList                                                                 ║
╚═════════════════════════════════════════════════════════════════════════════════════╝

╔═══════════════════════════════════════════════════════════╗
║ MERGESORT(List)                                           ║
║   if List consists of a single element                    ║
║     return List                                           ║
║   FirstHalf <- first half of List                         ║
║   SecondHalf <- second half of List                       ║
║   SortedFirstHalf <- MERGESORT(FirstHalf)                 ║
║   SortedSecondHalf <- MERGESORT(SecondHalf)               ║
║   SortedList <- MERGE(SortedFirstHalf, SortedSecondHalf)  ║
║   return SortedList                                       ║
╚═══════════════════════════════════════════════════════════╝
- k sorted array
  https://www.youtube.com/watch?v=yQ84lk-EXTQ
"""
#!/usr/bin/env python
import time


def merge(left, right):
    result = []
    while len(left) > 0 and len(right) > 0:
        if left[0] <= right[0]:
            result.append(left.pop(0))
        else:
            result.append(right.pop(0))
    while len(left) > 0:
        result.append(left.pop(0))
    while len(right) > 0:
        result.append(right.pop(0))
    return result


def mergesort(list):
    if len(list) == 1:
        return list
    left = list[0: len(list) // 2]
    right = list[len(list) // 2:]
    left = mergesort(left)
    right = mergesort(right)
    return  merge(left, right)


def abs(num):
    """
    int -> int
    return abs(num)
    """
    return num if num >=0 else -num


def greedy_sort(perm):
    """
    [a] -> [[a]]
    returns sequence of permutations
    TODO: improve this function
    """
    result = []
    comp = sorted(list(map(lambda x: x if x >= 0 else -x, perm)))  # list for comparison
    i = 0
    length = len(perm)
    while i < length:
        if abs(comp[i]) != perm[i]:
            idx = perm.index(comp[i]) if perm.count(comp[i]) else perm.index(-comp[i])
            perm = perm[:i] + list(map(lambda x: -x, perm[i:idx+1]))[::-1] + perm[idx+1:]
            result.append(perm)
            if perm[i] < 0:
                perm = perm[:i] + [-perm[i]] + perm[i+1:]
                result.append(perm)
        i += 1
    return result


def pretty_print(lls):
    for ls in lls:
        line = []
        for num in ls:
            if num >= 0:
                line.append('+' + str(num))
            else:
                line.append(str(num))
        print('(' + ' '.join(line) + ')')


def main():
    f = open('/home/wsl/rosalind/data/ba06a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    perm_str = lines[0].strip(')').strip('(').split()
    perm_int = list(map(int,lines[0].strip(')').strip('(').split()))
    f.close()

    start_time = time.time()
    lls = greedy_sort(perm_int)
    pretty_print(lls)

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
