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

solution by egeulgen
https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA4L.py
'''

import sys
MASS_TABLE = {'A': 71, 'C': 103, 'E': 129, 'D': 115, 'G': 57, 'F': 147, 'I': 113, 'H': 137, 'K': 128, 'M': 131,
              'L': 113, 'N': 114, 'Q': 128, 'P': 97, 'S': 87, 'R': 156, 'T': 101, 'W': 186, 'V': 99, 'Y': 163}


def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        temp = PrefixMass[i] + MASS_TABLE[Peptide[i]]
        PrefixMass.append(temp)
    LinearSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum


def linear_score(peptide, spectrum):
    pep_spec = LinearSpectrum(peptide)
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
        scores[i] = linear_score(peptide, spectrum)
    sorted_scores = sorted(scores.values(), reverse=True)
    threshold = sorted_scores[N - 1]
    return [leaderboard[idx] for idx, score in scores.items() if score >= threshold]


if __name__ == "__main__":
    '''
    Given: A leaderboard of linear peptides Leaderboard, a linear spectrum Spectrum, and an integer N.
    Return: The top N peptides from Leaderboard scored against Spectrum. Remember to use LinearScore.
    '''
    try:
        with open('/home/wsl/rosalind/data/ba04l.txt', 'r') as f:
            input_lines = [line.strip() for line in f.readlines()]
            Leaderboard = [line.strip() for line in input_lines[0].split(' ')]
            Spectrum = list(map(int, input_lines[1].split(' ')))
            N = int(input_lines[2].strip())
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    print(" ".join(Trim(Leaderboard, Spectrum, N)))
