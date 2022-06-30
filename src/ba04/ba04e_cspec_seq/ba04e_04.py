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
find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.

Return: Every amino acid string Peptide such that
        Cyclospectrum(Peptide) = Spectrum (if such a string exists).

Sample Dataset
0 113 128 186 241 299 314 427

Sample Output
186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

-------------------------------------------------

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


def Consistent(Peptide, Spectrum):
    if sum(Peptide) > Spectrum[-1] - MASSES[0]:
        return False
    spec = LinearSpectrum(Peptide)
    for mass in spec:
        if mass not in Spectrum:
            return False
    return True


def cyclopeptide_sequencing(spectrum):
    result = set()
    peptides = [[]]
    while peptides:
        peptides = expand(peptides)
        for peptide in peptides:
            if sum(peptide) == spectrum[-1]:
                if cyclospectrum_mass_peptide(peptide) == spectrum:
                    result.add("-".join(map(str, peptide)))
                peptides = [pep for pep in peptides if pep != peptide]
            elif not Consistent(peptide, spectrum):
                peptides = [pep for pep in peptides if pep != peptide]
    return result


if __name__ == "__main__":
    '''
    Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
    Return: Every amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a string exists).
    '''
    # spectrum = [int(x) for x in sys.stdin.readline().strip().split()]

    f = open('/home/wsl/rosalind/data/ba04e.txt', 'r')
    spectrum = [int(elem.strip()) for elem in f.readline().split()]
    f.close()

    print(" ".join(cyclopeptide_sequencing(spectrum)))