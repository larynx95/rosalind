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

TODO: Analyze this later.

solution by egeulgen
https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA4E.py
'''

import sys
MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def cyclospectrum_mass_peptide(peptide):
    spec = [0, sum(peptide)]
    temp = peptide + peptide
    for k in range(1, len(peptide)):
        for i in range(len(peptide)):
            subpeptide = temp[i:i + k]
            spec.append(sum(subpeptide))
    spec.sort()
    return spec


def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        temp = PrefixMass[i] + Peptide[i]
        PrefixMass.append(temp)
    LinearSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum


def expand(peptides):
    new_peptides = []
    for pep in peptides:
        for mass in MASSES:
            new_peptides.append(pep + [mass])
    return new_peptides


def Score(peptide, spectrum):
    pep_spec = cyclospectrum_mass_peptide(peptide)
    result = 0
    unique_masses = set(pep_spec + spectrum)
    for mass in unique_masses:
        result += min(pep_spec.count(mass), spectrum.count(mass))
    return result


def Trim(leaderboard, spectrum, N):
    if len(leaderboard) <= N:
        return leaderboard

    scores = {}
    for i, peptide in enumerate(leaderboard):
        scores[i] = Score(peptide, spectrum)

    sorted_scores = sorted(scores.values(), reverse=True)
    threshold = sorted_scores[N - 1]

    return [leaderboard[idx] for idx, score in scores.items() if score >= threshold]


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    leaderboard = [[]]
    leader_peptide = []

    while leaderboard:
        leaderboard = expand(leaderboard)
        for peptide in leaderboard:
            if sum(peptide) == spectrum[-1]:
                if Score(peptide, spectrum) > Score(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif sum(peptide) > spectrum[-1]:
                leaderboard = [pep for pep in leaderboard if pep != peptide]
        leaderboard = Trim(leaderboard, spectrum, N)
    return leader_peptide


if __name__ == "__main__":
    '''
    Given: An integer N and a collection of integers Spectrum.
    Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
    '''
    input_lines = sys.stdin.read().splitlines()
    N = int(input_lines[0])
    Spectrum = [int(x) for x in input_lines[1].split()]

    print("-".join(map(str, leaderboard_cyclopeptide_sequencing(Spectrum, N))))