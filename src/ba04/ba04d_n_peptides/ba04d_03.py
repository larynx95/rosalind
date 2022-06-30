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

═════════════════════════════════════════════════

https://github.com/Leoberium/BA/blob/master/Chapter4/BA4D.py
'''

aa_int_masses = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97,
    'V': 99, 'T': 101, 'C': 103, 'I': 113,
    'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137,
    'F': 147, 'R': 156, 'Y': 163, 'W': 186
}


def number_of_peptides(mass):
    arr = [0] * (mass + 1)  # index as peptide mass
    masses = set(aa_int_masses.values())
    # initialization
    for m in masses:
        arr[m] += 1
    # calculating the numbers
    min_mass = min(masses)
    for m in range(min_mass, mass + 1):
        for am in masses:
            r = m - am
            if r < 0:
                continue
            arr[m] += arr[r]
    return arr[mass]


def main():
    m = int(input())
    print(number_of_peptides(mass=m))


if __name__ == '__main__':
    main()