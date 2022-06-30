'''
Rosalind: BA4L
Trim a Peptide Leaderboard

The Trim algorithm, shown below,
sorts all peptides in Leaderboard according to their scores, resulting in a sorted Leaderboard.
Trim> then retains the "top N scoring peptides including ties",
and removes all other peptides from Leaderboard.

    Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
        for j <- 1 to |Leaderboard|
            Peptide <- j-th peptide in Leaderboard
            LinearScores(j) <- LinearScore(Peptide, Spectrum)
        sort Leaderboard according to the decreasing order of scores in LinearScores
        sort LinearScores in decreasing order
        for j <- N + 1 to |Leaderboard|
            if LinearScores(j) < LinearScores(N)
                remove all peptides starting from the j-th peptide from Leaderboard
            return Leaderboard
        return Leaderboard

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

solution by hadaarjan
'''

from collections import Counter

aminoacidMass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}


def linearSpectrum(peptide):
    """Input: An amino acid string Peptide.
     Output: The linear spectrum of Peptide."""
    prefixMass = [0]*((len(peptide)+1))
    for i in range(len(peptide)):
        prefixMass[i+1] = prefixMass[i] + aminoacidMass[peptide[i]]
    #print 'prefixMass', prefixMass
    linear_spectrum = [0]
    for i in range(len(prefixMass)-1):
        for j in range(i+1, len(prefixMass)):
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
    return [i[1] for i in sorted_scores if i[0] >= sorted_scores[int(N)-1][0]]


#Reading the iput from .txt file
f = open('/home/wsl/rosalind/data/ba04l.txt', 'r')
leaderboard, spectrum, n = [line.split() for line in f]
spectrum = map(int, spectrum)
N = int(n.pop())


#Calling the function
ans = trim_leaderboard(leaderboard, spectrum, N)


#Print the output in desired format
for i in ans:
    print(i)