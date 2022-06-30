"""
Rosalind: BA4L
Trim a Peptide Leaderboard

The Trim algorithm, shown below,
sorts all peptides in Leaderboard according to their scores, resulting in a sorted Leaderboard.
Trim> then retains the "top N scoring peptides including ties",
and removes all other peptides from Leaderboard.

╔══════════════════════════════════════════════════════════════════════════════════╗
║ Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)                         ║
║     for j <- 1 to |Leaderboard|                                                  ║
║         Peptide <- j-th peptide in Leaderboard                                   ║
║         LinearScores(j) <- LinearScore(Peptide, Spectrum)                        ║
║     sort Leaderboard according to the decreasing order of scores in LinearScores ║
║     sort LinearScores in decreasing order                                        ║
║     for j <- N + 1 to |Leaderboard|                                              ║
║         if LinearScores(j) < LinearScores(N)                                     ║
║             remove all peptides starting from the j-th peptide from Leaderboard  ║
║         return Leaderboard                                                       ║
║     return Leaderboard                                                           ║
╚══════════════════════════════════════════════════════════════════════════════════╝

Trim Problem
Trim a leaderboard of peptides.

Given: A leaderboard of linear peptides Leaderboard, a linear spectrum Spectrum, and an integer N.

Return: The top N peptides from Leaderboard scored against Spectrum. Remember to use LinearScore.

Sample Dataset
LAST ALST TLLT TQAS
0 71 87 101 113 158 184 188 259 271 372
2

Sample Output
LAST ALST

═════════════════════════════════════════════════

    [ Where am I? ]

    * Central Dogma
      -> translation (BA4A)
      -> reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      totally different from the Central Dogma
      ↓
    * Linear
      -> LinearSpectrum (BA4J)
      -> LinearScore (BA4K)
      -> the number of Linear peptides (BA4D)
      ↓
    * Cyclo
      -> Cyclospectrum (BA4C)
      -> CycloScore (BA4F)
      -> PREV: CycleoPeptideSequencing (BA4E)
      ↓
    * Leaderboard
      -> NEXT: LeaderboardCyclopeptideSequencing (BA4G)
        -> HERE: Trim (BA4L)

Plan 1.
- What's the meaning of the phrase, "top N scoring peptides including ties"?

  sample list: [(10,'a'),(8,'b'),(8,'z'),(3,'c'),(0,'d')]
  top 3 including ties    : [(10,'a'),(8,'b'),(8,'z'),(3,'c')]

  top 1: (10,'a')
  top 2: (8,'b'), (8,'z')
  top 3: (3,'c')

- example
  linear_score('LAST',[0,71,87,101,113,158,184,188,259,271,372])  --> 11
  linear_score('ALST',[0,71,87,101,113,158,184,188,259,271,372])  --> 9
  linear_score('TLLT',[0,71,87,101,113,158,184,188,259,271,372])  --> 5
  linear_score('TQAS',[0,71,87,101,113,158,184,188,259,271,372])  --> 5

  leaderboard = ['LAST','ALST','TLLT','TQAS']
  lspectrum   = [0,71,87,101,113,158,184,188,259,271,372]

  top 1: trim(leaderboard, lspectrum, 1)  --> (11,'LAST')
  top 2: trim(leaderboard, lspectrum, 2)  --> (9,'ALST')
  top 3: trim(leaderboard, lspectrum, 3)  --> (5,'TLLT'), (5,'TQAS')
  top 4: trim(leaderboard, lspectrum, 4)  --> ?

═════════════════════════════════════════════════

References:
- sorting a list of tuples
  How to sort a list/tuple of lists/tuples by the element at a given index?
  https://stackoverflow.com/questions/3121979/how-to-sort-a-list-tuple-of-lists-tuples-by-the-element-at-a-given-index
  sorted_by_second = sorted(data, key=lambda tup: tup[1])
- Python list sort in descending order
  https://stackoverflow.com/questions/4183506/python-list-sort-in-descending-order
  sorted(timestamps, reverse=True)
"""

#!/usr/bin/env python3
import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
    'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES = [57, 71, 87, 97, 99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]


def linear_spectrum(peptide):
    """
    string -> [int]
    algorithm in textbook
    >>> linear_spectrum_02('NQEL')
        [0,113,114,128,129,242,242,257,370,371,484]
    """
    # constructing cummulative PrefixMass list
    prefix_mass = [0]
    for i in range(len(peptide)):
        for j in range(len(AACIDS)):
            if peptide[i] == AACIDS[j]:
                prefix_mass.append(prefix_mass[-1] + AACID_MASSES[j])
    # constructing spectrum
    linear_spectrum = [0]
    for i in range(len(peptide)):
        for j in range(0, len(peptide)-i):
            linear_spectrum.append(prefix_mass[j+i+1] - prefix_mass[j])
    return sorted(linear_spectrum)


def linear_score(peptide, lspectrum):
    """
    string -> [int]
    >>> linear_score('NQEL', [0,99,113,114,128,227,257,299,355,356,370,371,484])
        8
    """
    common = []
    theoretical_lspectrum = linear_spectrum(peptide)
    for mass in lspectrum:
        if mass in theoretical_lspectrum:
            theoretical_lspectrum.remove(mass)
            common.append(mass)
    return len(common)


def trim(leader_board, spectrum, n):
    """
    >>> trim(['LAST','ALST','TLLT','TQAS'],[0,71,87,101,113,158,184,188,259,271,372],2)
        [(11, 'LAST'), (9, 'ALST')]
    """
    # get a sorted list of tuples (score, peptide)
    scores = []
    for peptide in leader_board:
        score = linear_score(peptide, spectrum)
        scores.append((score, peptide))
    scores.sort(reverse=True)
    print(scores)
    # trim
    for j in range(n, len(scores)):
        if scores[j][0] <= scores[n-1][0]:
            scores = scores[:j]
            return scores


def pretty_print(lts):
    for tup in lts:
        print(tup[1], end=' ')
    print()


def main():
    try:
        with open('/home/wsl/rosalind/data/ba04l.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            lb = [line.strip() for line in lines[0].split(' ')]
            lspec = list(map(int, lines[1].split(' ')))
            num = int(lines[2].strip())
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    result = trim(lb, lspec, num)
    pretty_print(result)
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()
