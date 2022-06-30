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

solution by hadaarjan
'''

aminoacid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
aminoacidMass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'L':113, 'N':114, 'D':115, 'K':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
def expand(leaderboard):
    """Expands each peptide/aminoacid in leaderboard by all 18 aminoacids with distinct masses."""
    expanded = []
    for i in leaderboard:
        expanded += [i+j for j in aminoacidMass.keys()]
    return expanded


def mass(peptide):
    """Calculates the mass of peptide using the aminoacidMass dictionary"""
    massOfPeptide = 0
    for i in peptide:
        massOfPeptide += aminoacidMass[i]
    return massOfPeptide


def cyclicSpectrum(peptide):
    """Input: An amino acid string Peptide.
     Output: The cyclic spectrum of Peptide."""
    prefixMass = [0]*((len(peptide)+1))
    for i in range(len(peptide)):
        prefixMass[i+1] = prefixMass[i] + aminoacidMass[peptide[i]]
    peptideMass = prefixMass[len(peptide)]
    cyclic_spectrum = [0]
    for i in range(len(prefixMass)-1):
        for j in range(i+1, len(prefixMass)):
            cyclic_spectrum.append(prefixMass[j] - prefixMass[i])
            if i > 0 and j < (len(prefixMass)-1):
                cyclic_spectrum.append(peptideMass - (prefixMass[j] - prefixMass[i]))
    return sorted(cyclic_spectrum)


from collections import Counter
def score_peptide(peptide, spectrum):
    """Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum)."""
    spectrum_peptide = cyclicSpectrum(peptide)
    c1, c2 = Counter(spectrum_peptide), Counter(spectrum)
    return sum([min(n, c2[k]) for k,n in c1.items()])


def linearSpectrum(peptide):
    """Input: An amino acid string Peptide.
     Output: The linear spectrum of Peptide."""
    prefixMass = [0]*((len(peptide)+1))
    for i in xrange(len(peptide)):
        prefixMass[i+1] = prefixMass[i] + aminoacidMass[peptide[i]]
    #print 'prefixMass', prefixMass
    linear_spectrum = [0]
    for i in xrange(len(prefixMass)-1):
        for j in xrange(i+1, len(prefixMass)):
            linear_spectrum.append(prefixMass[j] - prefixMass[i])
    return sorted(linear_spectrum)


def score_linear_peptide(peptide, spectrum):
    """Compute the score of a linear peptide with respect to a spectrum.
     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The linear score of Peptide with respect to Spectrum, LinearScore(Peptide, Spectrum)."""
    spectrum_linear_peptide = linearSpectrum(peptide)
    c3, c4 = Counter(spectrum_linear_peptide), Counter(spectrum)
    return sum([min(n, c4[k]) for k,n in c3.items()])


def trim_leaderboard(leaderboard, spectrum, N):
    """Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
     Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum."""
    scores =  [[score_linear_peptide(peptide, spectrum), peptide] for peptide in leaderboard]
    sorted_scores = sorted(scores, reverse = True)
    if len(leaderboard) <= N:
        return [i[1] for i in sorted_scores]
    else:
        return [i[1] for i in sorted_scores if i[0] >= sorted_scores[int(N)-1][0]]


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    """ Input: An integer N and a collection of integers Spectrum.
     Output: LeaderPeptide after running LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)"""
    leaderboard = aminoacid
    leaderpeptide = ''
    parentmass = max(spectrum)
    while len(leaderboard) > 0:
        leaderboard = expand(leaderboard)
        for peptide in leaderboard[:]:
            if mass(peptide) == parentmass:
                if score_peptide(peptide, spectrum) > score_peptide(leaderpeptide, spectrum):
                    leaderpeptide = peptide
            elif mass(peptide) > parentmass:
                leaderboard.remove(peptide)
        leaderboard = trim_leaderboard(leaderboard, spectrum, N)
    return leaderpeptide