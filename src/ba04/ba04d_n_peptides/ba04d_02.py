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

https://github.com/andreyrozumnyi/rosalind/blob/master/Chapter%204/BA4D.py
'''

""" Compute the Number of Peptides of Given Total Mass """

keys = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
def calc_varinats(m):
    ways = [0]*(m + 1)
    index = m
    ways[m] = 1
    while index > 0:
        for key in keys:
            ways[index-key] += ways[index]
        index -= 1
        while ways[index] == 0:
            index -= 1
    return ways[0]

if __name__ == "__main__":
    print(calc_varinats(1024))