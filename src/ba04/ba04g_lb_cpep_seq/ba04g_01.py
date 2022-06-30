'''
Rosalind: BA4G
Implement LeaderboardCyclopeptideSequencing

We have thus far worked with theoretical spectra of cyclic peptides,
in which the mass of every subpeptide is given.
This inflexibility presents a practical barrier,
since mass spectrometers generate spectra that are far from ideal
— they are characterized by having both false masses and missing masses.
A false mass is present in the experimental spectrum but absent from the theoretical spectrum;
a missing mass is present in the theoretical spectrum but absent from the experimental spectrum (see Figure 1).

To generalize "Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum"
to handle "noisy" spectra having false and missing masses,
we need to relax the requirement that a candidate peptide’s theoretical spectrum must match the experimental spectrum exactly,
and instead incorporate a scoring function that will select the peptide
whose theoretical spectrum matches the given experimental spectrum the most closely.
Given a cyclic peptide Peptide and a spectrum Spectrum,
we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum.
Recalling Figure 1, if

Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},
then Score("NQEL", Spectrum) = 11.

To limit the number of candidate peptides under consideration,
we will use a Leaderboard, which holds the N highest scoring candidate peptides for further extension.
At each step, we will expand all candidate peptides found in Leaderboard by adding every possible amino acid to the end.
Then, we will eliminate those peptides whose newly calculated scores are not high enough to keep them on the Leaderboard.
This idea is similar to the notion of a “cut” in a golf tournament;
after the cut, only the top N golfers are allowed to play in the next round,
since they are the only players who have a reasonable chance of winning.

To be fair, a cut should include anyone who is tied with the Nth-place competitor.
Thus, Leaderboard should be trimmed down to the “N highest-scoring peptides including ties”,
which may include more than N peptides.
Given a list of peptides Leaderboard, a spectrum Spectrum, and an integer N, Cut(Leaderboard, Spectrum, N)
returns the top N highest-scoring peptides in <Leaderboard (including ties) with respect to Spectrum.
We now introduce LEADERBOARDCYCLOPEPTIDESEQUENCING.
In what follows, the 0-peptide is the peptide "" containing no amino acids.

    LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
        Leaderboard <- {0-peptide}
        LeaderPeptide <- 0-peptide
        while Leaderboard is non-empty
            Leaderboard <- Expand(Leaderboard)
            for each Peptide in Leaderboard
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                        LeaderPeptide <- Peptide
                else if Mass(Peptide) > ParentMass(Spectrum)
                    remove Peptide from Leaderboard
            Leaderboard <- Cut(Leaderboard, Spectrum, N)
        output LeaderPeptide

Implement LeaderboardCyclopeptideSequencing
Given: An integer N and a collection of integers Spectrum.

Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).

Sample Dataset
10
0 71 113 129 147 200 218 260 313 331 347 389 460

Sample Output
113-147-71-129

═════════════════════════════════════════════════

solution by Leoberium
https://github.com/Leoberium/BA/blob/master/Chapter4/BA4G.py
'''

#!/usr/bin/env python3
import sys


masses = [
    57, 71, 87, 97, 99,
    101, 103, 113, 114, 115,
    128, 129, 131, 137, 147,
    156, 163, 186
]


def mass(peptide):
    m = sum(map(int, peptide.split('-')))
    return m


def expand(ps):
    for peptide in ps.copy():
        ps.remove(peptide)
        for m in masses:
            if peptide:
                new_peptide = peptide + '-' + str(m)
            else:
                new_peptide = str(m)
            ps.add(new_peptide)
    return ps


def linear_spectrum(peptide):
    cs = [0]
    if not peptide:
        return cs
    composition = list(map(int, peptide.split('-')))
    n = len(composition)
    for k in range(1, n):
        for i in range(n - k + 1):
            cs.append(sum(composition[i:i+k]))
    cs.append(mass(peptide))
    cs.sort()
    return cs


# linear score to trim the leaderboard
def linear_score(peptide, spectrum):
    peptide_spectrum = linear_spectrum(peptide)
    p_dict, s_dict = dict(), dict()
    for value in peptide_spectrum:
        if value in p_dict:
            p_dict[value] += 1
        else:
            p_dict[value] = 1
    for value in spectrum:
        if value in s_dict:
            s_dict[value] += 1
        else:
            s_dict[value] = 1
    s = 0
    for value in p_dict:
        if value in s_dict:
            s += min(p_dict[value], s_dict[value])
    return s


def cyclic_spectrum(peptide):
    cs = [0]
    if not peptide:
        return cs
    composition = list(map(int, peptide.split('-')))
    n = len(composition)
    for k in range(1, n):
        for i in range(n):
            if i + k <= n:
                cs.append(sum(composition[i:i+k]))
            else:
                r = i + k - n
                cs.append(sum(composition[i:] + composition[:r]))
    cs.append(mass(peptide))
    cs.sort()
    return cs


# cyclic score for the leader peptide
def cyclic_score(peptide, spectrum):
    peptide_spectrum = cyclic_spectrum(peptide)
    p_dict, s_dict = dict(), dict()
    for value in peptide_spectrum:
        if value in p_dict:
            p_dict[value] += 1
        else:
            p_dict[value] = 1
    for value in spectrum:
        if value in s_dict:
            s_dict[value] += 1
        else:
            s_dict[value] = 1
    s = 0
    for value in p_dict:
        if value in s_dict:
            s += min(p_dict[value], s_dict[value])
    return s


def trim(leader_board, spectrum, n):
    scores = []
    for peptide in leader_board:
        score = linear_score(peptide, spectrum)
        scores.append((score, peptide))
    scores.sort(reverse=True)
    for j in range(n, len(scores)):
        if scores[j][0] < scores[n-1][0]:
            scores = scores[:j]
            break
    leader_board = set()
    for score, peptide in scores:
        leader_board.add(peptide)
    return leader_board


def leaderboard_cyclopeptide_sequencing(spectrum, n):
    leader_board = {''}
    leader_peptide = ''
    parent_mass = max(spectrum)
    while leader_board:
        leader_board = expand(leader_board)
        for peptide in leader_board.copy():
            if mass(peptide) == parent_mass:
                if cyclic_score(peptide, spectrum) > \
                        cyclic_score(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif mass(peptide) > parent_mass:
                leader_board.remove(peptide)
        leader_board = trim(leader_board, spectrum, n)
    return leader_peptide


def main():
    n = int(sys.stdin.readline())
    spectrum = list(map(int, sys.stdin.readline().split()))
    print(leaderboard_cyclopeptide_sequencing(spectrum, n))


if __name__ == '__main__':
    main()