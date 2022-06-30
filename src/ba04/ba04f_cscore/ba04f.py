"""
Rosalind: BA4F
Compute the Score of a Cyclic Peptide Against a Spectrum

To generalize the Cyclopeptide Sequencing Problem
from "Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum" to handle noisy spectra,
we need to relax the requirement
that a candidate peptide's theoretical spectrum must match the experimental spectrum exactly,
and instead incorporate a scoring function that will select the peptide
whose theoretical spectrum matches the given experimental spectrum the most closely.
Given a cyclic peptide Peptide and a spectrum Spectrum,
we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum.
Recalling our example above, if

> Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484}, then Score(NQEL, Spectrum) = 11.

The scoring function should take into account the multiplicities of shared masses,
i.e., how many times they occur in each spectrum. For example,
suppose that Spectrum is the theoretical spectrum of NQEL;
for this spectrum, mass 242 has multiplicity 2.
If 242 has multiplicity 1 in the theoretical spectrum of Peptide,
then 242 contributes 1 to Score(Peptide, Spectrum).
If 242 has larger multiplicity in the theoretical spectrum of Peptide,
then 242 contributes 2 to Score(Peptide, Spectrum).

Cyclic Peptide Scoring Problem
Compute the score of a cyclic peptide against a spectrum.

Given: An amino acid string Peptide and a collection of integers Spectrum.

Return: The score of Peptide against Spectrum, Score(Peptide, Spectrum).

Sample Dataset
NQEL
0 99 113 114 128 227 257 299 355 356 370 371 484

Sample Output
11

═════════════════════════════════════════════════

    [ Where am I? ] How do we sequence antibiotics?

    * Central Dogma
      -> translation (BA4A)
      -> reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      -> Linear/Cyclo - spectrum (BA4C)
      -> the Number of Linear Peptides of Given Total Mass (BA4D)
      -> PREV: CycloPeptideSequencing (BA4E)
      -> HERE: Score(Peptide, Spectrum) (BA4F)

Plan 1.
       L   N   Q   E   LN  NQ  EL  QE      LNQ ELN QEL NQE NQEL
  0    113 114 128 129 227 242 242 257     355 356 370 371 484  <-- cyclic spectrum of NQEL
  0 99 113 114 128     227         257 299 355 356 370 371 484  <-- given spectrum
  1    2   3   4       5           6       7   8   9   10  11   <-- score

═════════════════════════════════════════════════

    [ Where am I? ]

    * Central Dogma
      -> translation (BA4A)
      -> reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      totally different from the Central Dogma
      ↓
    * Linear
      -> LinearSpectrum (BA4J)
      -> LinearScore (BA4K)
      -> the number of Linear peptides (BA4D)
      ↓
    * Cyclo
      -> PREV: Cyclospectrum (BA4C)
      -> HERE: CycloScore (BA4F)
      -> NEXT: CycleoPeptideSearch (BA4E)

References:
- Intersection of two lists including duplicates?
  https://stackoverflow.com/questions/37645053/intersection-of-two-lists-including-duplicates
- Intersection of two lists, keeping duplicates in the first list
  https://stackoverflow.com/questions/26663371/intersection-of-two-lists-keeping-duplicates-in-the-first-list
"""

import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
    'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES = [57, 71, 87, 97, 99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]


def cyclo_spectrum(peptide):
    """
    algorithm in textbook
    >>> cyclo_spectrum('NQEL')
        [0,113,114,128,129,227,242,242,257,355,356,370,371,484]
    """
    prefix_mass = [0]
    for aacid in peptide:
        prefix_mass.append(prefix_mass[-1] + WT_TABLE[aacid])
    cspectrum = [0]
    peptide_mass = prefix_mass[-1]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            cspectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cspectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cspectrum)


def cyclo_spectrum_ints(ipeptide):
    """
    [int] -> [int]
    returns a cyclo-spectrum from a list of single amino acid masses
    >>> cyclo_spectrum_ints([113,128,186])
        [0,113,128,186,241,299,314,427]
    """
    prefix_mass = [0]
    for mass in ipeptide:
        prefix_mass.append(prefix_mass[-1] + mass)
    cspectrum = [0]
    peptide_mass = prefix_mass[-1]
    for i in range(len(ipeptide)):
        for j in range(i+1, len(ipeptide)+1):
            cspectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(ipeptide):
                cspectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cspectrum)


def cyclo_score(peptide, spectrum):
    """
    >>> cyclo_score('NQEL', [0,99,113,114,128,227,257,299,355,356,370,371,484])
        11
    """
    common = []
    theoretical_cspectrum = cyclo_spectrum(peptide)
    for mass in theoretical_cspectrum:
        if mass in spectrum:
            spectrum.remove(mass)
            common.append(mass)
    return len(common)


def cyclo_score_ints(ipeptide, cspectrum):
    """
    >>> cyclo_score_ints([114-128-129-113], [0,99,113,114,128,227,257,299,355,356,370,371,484])
    """
    common = []
    theoretical_cspectrum = cyclo_spectrum_ints(ipeptide)
    for mass in cspectrum:
        if mass in theoretical_cspectrum:
            theoretical_cspectrum.remove(mass)
            common.append(mass)
    return len(common)


def main():
    f = open('/home/wsl/rosalind/data/ba04f.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    peptide = lines[0]
    cspectrum = [int(line.strip()) for line in lines[1].split(' ')]
    f.close()

    start_time = time.time()
    sc = cyclo_score(peptide, cspectrum)
    print(sc)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
