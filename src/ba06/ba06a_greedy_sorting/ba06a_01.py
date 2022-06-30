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

Info.
  * identity permutation:
    ; ordered from smallest to largest with positive directions

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

═════════════════════════════════════════════════

References:
- solution by Matthew Burch (cryptc)

"""

import time


def PrintP(P):
    print("({})".format(" ".join(map(lambda k: "{0:+}".format(k),P))))


def StitchTogether(P, k, idx):
    return P[:k] + Reversal(P[k:idx+1]) + P[idx+1:]


def FindIndexFromValue(P, val):
    if val in P:
        return P.index(val)
    return P.index(-val)


def Reversal(arr):
    return [-k for k in reversed(arr)]


def GreedySorting(signed_permutation_P):
    P = [int(i) for i in signed_permutation_P[1:-1].split()]
    for k in range(len(P)):
        idx = FindIndexFromValue(P, k+1)
        if P[k] != abs(P[k]) or k != idx:
            P = StitchTogether(P, k, idx)
            PrintP(P)
            if P[k] != abs(P[k]):
                P[k] = abs(P[k])
                PrintP(P)


def main():
    f = open('/home/wsl/rosalind/data/ba06a.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    perm_str = lines[0].strip(')').strip('(').split()
    perm_int = list(map(int,lines[0].strip(')').strip('(').split()))
    f.close()

    start_time = time.time()
    GreedySorting(lines[0])
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
