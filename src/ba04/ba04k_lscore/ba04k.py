"""
Rosalind: BA4K
Compute the Score of a Linear Peptide

Linear Peptide Scoring Problem
Compute the score of a linear peptide with respect to a spectrum.

Given: An amino acid string Peptide and a collection of integers LinearSpectrum.

Return: The linear score of Peptide against Spectrum, LinearScore(Peptide, Spectrum).

Sample Dataset
NQEL
0 99 113 114 128 227 257 299 355 356 370 371 484

Sample Output
8

═════════════════════════════════════════════════

    [ Where am I? ] How do we sequence antibiotics?

    * Central Dogma
      -> translation (BA4A)
      -> reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      totally different from the Central Dogma
      ↓
    * Linear
      -> PREV: LinearSpectrum (BA4J)
      -> HERE: LinearScore (BA4K)
      -> NEXT: the number of Linear peptides (BA4D)

Plan 1.
- two lists

  theoretical : [0    113 114 128 129     242 242 257             370 371 484]
  experimental: [0 99 113 114 128     227         257 299 355 356 370 371 484]

  [0    113 114 128 129     242 242 257             370 371 484]
  [0 99 113 114 128     227         257 299 355 356 370 371 484]

- problems:
  a. There're duplicates in one or both lists!
     => element should be removed during comparison
  b. One of the lists can be modified during iteration!
     => select theoretical linear-spectrum as modifying one
     => TODO: I found weird thing. (below) Why???

═════════════════════════════════════════════════

References:
"""

import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
    'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES = [57, 71, 87, 97, 99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]


def linear_spectrum(peptide):
    """
    string -> [int]
    algorithm in textbook
    >>> linear_spectrum_02('NQEL')
        [0,113,114,128,129,242,242,257,370,371,484]
    """
    # constructing cummulative PrefixMass list
    prefix_mass = [0]
    for i in range(len(peptide)):
        for j in range(len(AACIDS)):
            if peptide[i] == AACIDS[j]:
                prefix_mass.append(prefix_mass[-1] + AACID_MASSES[j])
    # constructing spectrum
    linear_spectrum = [0]
    for i in range(len(peptide)):
        for j in range(0, len(peptide)-i):
            linear_spectrum.append(prefix_mass[j+i+1] - prefix_mass[j])
    return sorted(linear_spectrum)


def linear_score(peptide, lspectrum):
    """
    (string,[int]) -> int
    >>> linear_score('NQEL', [0,99,113,114,128,227,257,299,355,356,370,371,484])
        8
    """
    common = []
    theoretical_lspectrum = linear_spectrum(peptide)
    for mass in theoretical_lspectrum:
        if mass in lspectrum:
            lspectrum.remove(mass)  # <-- lspectrum is modified here
            common.append(mass)
    return len(common)


def linear_score_better(peptide, lspectrum):
    """
    (string,[int]) -> int
    better version, cus lspectrum is intact
    >>> linear_score_better('NQEL', [0,99,113,114,128,227,257,299,355,356,370,371,484])
        8
    """
    common = []
    theoretical_lspectrum = linear_spectrum(peptide)
    for mass in lspectrum:
        if mass in theoretical_lspectrum:
            theoretical_lspectrum.remove(mass)  # <-- theoretical_lspectrum is modifed here, instead of lspectrum
            common.append(mass)
    return len(common)


def linear_score_WEIRD(peptide, lspectrum):
    """
    (string,[int]) -> int
    I copied the list to prevent the original list from being changed during for-loop.
    However, strange results came out!!
    >>> linear_score_WEIRD('NQEL', [0,99,113,114,128,227,257,299,355,356,370,371,484])
        8
    Above result is OK, but this function returns weird result with downloaded dataset!
    TODO: Why? Fix this!
    """
    common = []
    theoretical_lspectrum = linear_spectrum(peptide)
    for mass in theoretical_lspectrum:
        copied = lspectrum.copy()  # <-- something mysterious thing happened here!
        if mass in copied:
            copied.remove(mass)
            common.append(mass)
    return len(common)


def main():
    try:
        with open('/home/wsl/rosalind/data/ba04k.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            peptide = lines[0]
            lspec = list(map(int, lines[1].split(' ')))
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    print(linear_score(peptide, lspec))
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()
