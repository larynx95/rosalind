'''
Rosalind: BA4D
Compute the Number of Peptides of Given Total Mass

    Counting Peptides with Given Mass Problem:
        Compute the number of peptides of given mass.
        Input: An integer m.
        Output: The number of "linear" peptides having integer mass m.
    Hint.
        Not easy! Use dynamic programming.

In "Generate the Theoretical Spectrum of a Cyclic Peptide",
we generated the theoretical spectrum of a known cyclic peptide.
Although this task is relatively easy, our aim in mass spectrometry is to solve the reverse problem:
we must reconstruct an unknown peptide from its experimental spectrum.
We will start by assuming that a biologist is lucky enough to generate an ideal experimental spectrum 'Spectrum',
which is one coinciding with the peptide’s theoretical spectrum.
Can we reconstruct a peptide whose theoretical spectrum is 'Spectrum'?

Denote the total mass of an amino acid string Peptide as Mass(Peptide).
In mass spectrometry experiments, whereas the peptide that generated a spectrum is unknown,
the peptide’s mass is typically known and is denoted "ParentMass(Spectrum)".
Of course, given an ideal experimental spectrum,
Mass(Peptide) is given by the largest mass in the spectrum.

A brute force approach to reconstructing a peptide from its theoretical spectrum
would generate all possible peptides whose mass is equal to "ParentMass(Spectrum)"
and then check which of these peptides has theoretical spectra matching Spectrum.
However, we should be concerned about the running time of such an approach:
how many peptides are there having mass equal to "ParentMass(Spectrum)"?

Counting Peptides with Given Mass Problem
Compute the number of peptides of given total mass.

Given: An integer m.

Return: The number of 'linear' peptides having integer mass m.

Sample Dataset
1024

Sample Output
14712706211

-------------------------------------------------

https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA4D.py
'''

import sys
MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

def CountPeptides(Mass):
    NumPeptides={}
    for i in range(57):
        NumPeptides[i] = 0
    for mass in range(57, Mass + 1):
        NumPeptides[mass] = MASSES.count(mass)
        for int_mass in MASSES:
            if mass >= int_mass:
                if NumPeptides[mass - int_mass] > 0:
                    NumPeptides[mass] += NumPeptides[mass - int_mass]
    return NumPeptides[Mass]


if __name__ == "__main__":
    '''
    Given: An integer m.
    Return: The number of linear peptides having integer mass m.
    '''
    m = int(sys.stdin.readline().strip())

    print(CountPeptides(m))