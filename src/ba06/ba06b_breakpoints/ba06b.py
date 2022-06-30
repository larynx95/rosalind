"""
Rosalind: BA6B
Compute the Number of Breakpoints in a Permutation

Number of Breakpoints Problem
Find the number of breakpoints in a permutation.

Given: A signed permutation P.

Return: The number of breakpoints in P.

Sample Dataset
(+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14)

Sample Output
8

═════════════════════════════════════════════════

Info.
  * Adjacencies(P), Breakpoints(P)
    (0, p1, ... , pn, (n+1))

    (+1 +2 ... +n)        : all adjacencies
    (-n -(n-1) ... -2 -1) : all adjacencies except (0 -n), (-1 (n+1))

  * sample dataset
      +3   +4   +5   -12   -8  -7  -6   +1   +2   +10   +9   -11   +13   +14
    | +3   +4   +5 | -12 | -8  -7  -6 | +1   +2 | +10 | +9 | -11 | +13   +14
       v
    (  0, +3)    breakpoint   + 1
    ( +3, +4)    adjacency
    ( +4, +5)    adjacency
    ( +5,-12)    breakpoint   + 2
    (-12, -8)    breakpoint   + 3
    ( -8, -7)    adjacency
    ( -7, -6)    adjacency
    ( -6, +1)    breakpoint   + 4
    ( +1, +2)    adjacency
    ( +2,+10)    breakpoint   + 5
    (+10, +9)    breakpoint   + 6
    ( +9,-11)    breakpoint   + 7
    (-11,+13)    breakpoint   + 8   the number of breakpoints is 8
    (+13,+14)    adjacency
    (+14,+15)    adjacency
          ^

  * example in textbook
                                                                              breakpoints
    0  |+3   +4   +5 | 12 |  8    7    6 | +1   +2 |+10 | +9 | 11 |+13  +14   8
    0  |+3   +4   +5 |+11 |  9 | 10 |  2    1 | +6   +7   +8 |+12  +13  +14   7
    0   +1   +2 |+10 | +9 | 11 |  5    4    3 | +6   +7   +8 |+12  +13  +14   6
    0   +1   +2   +3   +4   +5 |+11 |  9 | 10 | +6   +7   +8 |+12  +13  +14   5
    0   +1   +2   +3   +4   +5 | +9 | 11   10 | +6   +7   +8 |+12  +13  +14   4
    0   +1   +2   +3   +4   +5 | +9 |  8    7    6 |+10  +11  +12  +13  +14   3
    0   +1   +2   +3   +4   +5   +6   +7   +8  | 9 |+10  +11  +12  +13  +14   2
    0   +1   +2   +3   +4   +5   +6   +7   +8   +9  +10  +11  +12  +13  +14   0

Plan 1.
  * prepend 0 to the list, append n+1 to the list
  * p[n+1] - p[n] != 1, then ++break

═════════════════════════════════════════════════

References:
-
"""
#!/usr/bin/env python
import time


def num_breakpoints(chromosome):
    """
    [int] -> int
    returns the number of breakpoints
    >>> num_breakpoints([+3,+4,+5,-12,-8,-7,-6,+1,+2,+10,+9,-11,+13,+14])
        8
    """
    count = 0
    chromosome = [0] + chromosome + [len(chromosome)+1]  # len(P)+1, not len(P)
    for i in range(1, len(chromosome)):
        if chromosome[i] - chromosome[i-1] != 1:
            count += 1
    return count


def main():
    f = open('../data/ba6b.txt', 'r')
    line = f.readline().strip('\n').strip()
    perm_str = line[1:-1].split(' ')
    chromosome = list(map(int, perm_str))
    f.close()

    start_time = time.time()
    print(num_breakpoints(chromosome))

    # solution by other person - 'root' WOW! Just one line!
    print(sum(map(lambda x,y: x-y != 1, chromosome+[len(chromosome)+1], [0]+chromosome)))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()