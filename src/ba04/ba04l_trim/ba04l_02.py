"""
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

solution by Elmar Hinz
"""

def TrimLeaderboard(board, spectrum, size):
    map = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101,
           "C": 103, "I": 113, "L": 113, "N": 114, "D": 115, "K": 128,
           "Q": 128, "E": 129, "M": 131, "H": 137, "F": 147, "R": 156,
           "Y": 163, "W": 186}

    def linearSpectrum(peptide):
        spectrum = [0]
        for s in range(len(peptide)):
            for l in range(1, len(peptide) - s + 1):
                spectrum.append(sum((peptide)[s: s + l]))
        return sorted(spectrum)

    def score(peptide):
        score, spectrum2 = 0, linearSpectrum([map[c] for c in peptide])
        for mass in spectrum:
            if mass in spectrum2:
                spectrum2.remove(mass)
                score += 1
        return score

    def trim():
        scores = [score(peptide) for peptide in board]
        sorting = sorted(range(len(scores)),
                key=lambda k: scores[k], reverse = True)
        sortedBoard = [board[i] for i in sorting]
        sortedScores = [scores[i] for i in sorting]
        for i in range(size, len(board)):
            if sortedScores[i] < sortedScores[size - 1]:
                return sortedBoard[0:i]

    return trim()