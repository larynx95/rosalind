"""
Rosalind: BA4H
Generate the Convolution of a Spectrum

We define the convolution of a cyclic spectrum by taking all positive differences of masses in the spectrum.
Figure 1 shows the convolution of the theoretical spectrum of "NQEL".

As predicted, some of the values in Figure 2.12 appear more frequently than others.
For example, 113 (the mass of "L") appears eight times in the convolution of the theoretical spectrum of "NQEL";
we say that 113 has multiplicity 8.
Six of the eight occurrences of 113 correspond to subpeptide pairs
differing in an "L": "L" and ""; "LN" and "N"; "EL" and "E"; "LNQ" and "NQ"; "QEL" and "QE"; "NQEL" and "NQE".

Spectral Convolution Problem
Compute the convolution of a spectrum.

Given: A collection of integers Spectrum.

Return: The list of elements in the convolution of Spectrum in decreasing order of their multiplicities.
If an element has multiplicity k, it should appear exactly k times.

Sample Dataset
0 137 186 323

Sample Output
137 137 186 186 323 49

═════════════════════════════════════════════════

Plan 1.

References:
"""

import time
import sys


def convolution_wrong(spectrum):
    spectrum.sort()
    freq_dict = {}
    for i in range(1, len(spectrum)):
        for j in range(0, i):
            diff = spectrum[i] - spectrum[j]
            if diff not in freq_dict:  freq_dict[diff] = 1
            else: freq_dict[diff] += 1
    sorted_mass_list = [k for k, _ in sorted(freq_dict.items(), key=lambda item: item[1], reverse=True)]
    conv = []
    for mass in sorted_mass_list:
        conv += [mass] * freq_dict[mass]
    return conv


def convolution(spectrum):
    spectrum.sort()
    conv = []
    for i in range(len(spectrum) - 1):
        for j in range(i, len(spectrum)):
            if spectrum[j] - spectrum[i] != 0:
                conv.append(spectrum[j] - spectrum[i])
    freq_dict = {}
    for mass in set(conv):
        freq_dict[mass] = conv.count(mass)
    sorted_mass_list = [k for k, _ in sorted(freq_dict.items(), key=lambda item: item[1], reverse=True)]
    conv = []
    for mass in sorted_mass_list:
        conv += [mass] * freq_dict[mass]
    return conv


if __name__ == "__main__":
    f = open('/home/wsl/rosalind/data/ba04h.txt', 'r')
    spectrum = [int(x) for x in f.readline().strip().split()]
    f.close()
    print(" ".join(map(str, convolution(spectrum))))

    #with open('../data/ba4h_output.txt', 'w') as f:
    #    original_stdout = sys.stdout
    #    sys.stdout = f # Change the standard output to the file we created.
    #    print(" ".join(map(str, convolution(spectrum))))
    #    sys.stdout = original_stdout
