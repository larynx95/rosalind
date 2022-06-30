'''
Rosalind: BA4E
Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

In "Compute the Number of Peptides of Given Total Mass",
we first encountered the problem of reconstructing a cyclic peptide from its theoretical spectrum;
this problem is called the Cyclopeptide Sequencing Problem and is given below.
It is solved by the following algorithm.

    CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides <- a set containing only the empty peptide
        while Peptides is nonempty
            Peptides <- Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides

Cyclopeptide Sequencing Problem
Given an ideal experimental spectrum,
find a "cyclic" peptide whose theoretical spectrum matches the experimental spectrum.

Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.

Return: Every amino acid string Peptide such that
        Cyclospectrum(Peptide) = Spectrum (if such a string exists).

Sample Dataset
0 113 128 186 241 299 314 427

Sample Output
186-128-113  186-113-128  128-186-113  128-113-186  113-186-128  113-128-186

-------------------------------------------------

solution by Leoberium
https://github.com/Leoberium/BA/blob/master/Chapter4/BA4E.py
'''

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


def cyclospectrum(peptide):
    '''
    >>> cyclospectrum('57-113')
        [0,57,113,170]
    '''
    cs = [0]
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

print(cyclospectrum('57-113'))


def linear_spectrum(peptide):
    cs = [0]
    composition = list(map(int, peptide.split('-')))
    n = len(composition)
    for k in range(1, n):
        for i in range(n - k + 1):
            cs.append(sum(composition[i:i+k]))
    cs.append(mass(peptide))
    cs.sort()
    return cs


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


def not_consistent(peptide, spectrum):
    peptide_spectrum = linear_spectrum(peptide)
    for value in peptide_spectrum:
        if value not in spectrum:
            return True
    return False


def cyclopeptide_sequencing(spectrum):
    peptides = {''}
    matches = set()
    parent_mass = max(spectrum)
    while peptides:
        peptides = expand(peptides)
        for peptide in peptides.copy():
            if mass(peptide) == parent_mass:
                if cyclospectrum(peptide) == spectrum:
                    matches.add(peptide)
                peptides.remove(peptide)
            elif not_consistent(peptide, spectrum):
                peptides.remove(peptide)
    return matches


def main():
    # spectrum = list(map(int, sys.stdin.readline().split()))

    f = open('/home/wsl/rosalind/data/ba04e.txt', 'r')
    spectrum = [int(elem.strip()) for elem in f.readline().split()]
    f.close()

    print(*cyclopeptide_sequencing(spectrum))


if __name__ == '__main__':
    main()