"""
Rosalind: BA4J
Generate the Theoretical Spectrum of a Linear Peptide

Given an amino acid string Peptide, we will begin by assuming that it represents a linear peptide.
Our approach to generating its theoretical spectrum is based on the assumption
that the mass of any subpeptide is equal to the difference between the masses of two prefixes of Peptide.
We can compute an array PrefixMass storing the masses of each prefix of Peptide in increasing order, e.g.,
for Peptide = NQEL, PrefixMass = (0, 114, 242, 371, 484).
Then, the mass of the subpeptide of Peptide beginning at position i + 1
and ending at position j can be computed as PrefixMass(j) - PrefixMass(i).
For example, when Peptide = NQEL,

Mass(QE) = PrefixMass(3) - PrefixMass(1) = 371 - 114 = 257.
The pseudocode shown on the next step implements this idea.
It also represents the alphabet of 20 amino acids and their integer masses
as a pair of 20-element arrays AminoAcid and AminoAcidMass,
corresponding to the top and bottom rows of the following integer mass table, respectively.

Figure 1:
╔════════════════════════════════════════════════════════════════════════╗
║ LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)                      ║
║     PrefixMass(0) <- 0                                                 ║
║     for i <- 1 to |Peptide|                                            ║
║         for j <- 1 to 20                                               ║
║             if AminoAcid(j) =  i-th amino acid in Peptide              ║
║                 PrefixMass(i) <- PrefixMass(i - 1) + AminoAcidMass(j)  ║
║     LinearSpectrum <- a list consisting of the single integer 0        ║
║     for i <- 0 to |Peptide| - 1                                        ║
║         for j <- i + 1 to |Peptide|                                    ║
║             add PrefixMass(j) - PrefixMass(i) to LinearSpectrum        ║
║     return sorted list LinearSpectrum                                  ║
╚════════════════════════════════════════════════════════════════════════╝

Linear Spectrum Problem
Generate the ideal linear spectrum of a peptide.

Given: An amino acid string Peptide.

Return: The linear spectrum of Peptide.

Sample Dataset
NQEL

Sample Output
0 113 114 128 129 242 242 257 370 371 484

═════════════════════════════════════════════════

    [ Where am I? ] How do we sequence antibiotics?

    * Central Dogma
      -> translation (BA4A)
      -> PREV: reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      totally different from the Central Dogma
      ↓
    * Linear
      -> HERE: LinearSpectrum (BA4J)
      -> NEXT: LinearScore (BA4K)

Plan 1.
- using textbook algorithm
- Wouldn't it be better for this problem to be placed before 'ba4c'?

References:
"""

import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
       'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES = [57, 71, 87, 97, 99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]


def get_mass(peptide):
    """
    implementation of 'Mass(peptide)' function
    returns mass of a peptide
    """
    return sum([WT_TABLE[aacid] for aacid in peptide])


def linear_spectrum_01(peptide):
    """
    string -> [int]
    returns sorted linear spectrum of a peptide
    >>> linear_spectrum_01('NQEL')
        [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    """
    # get all subpeptides
    subpeptides = []
    for i in range(len(peptide)):
        for j in range(len(peptide) - i):
            subpeptides.append(peptide[j:i+j+1])
    # get spectrum from a list of all subpeptides
    spectrum = [0]
    for pep in subpeptides:
        spectrum.append(get_mass(pep))
    return sorted(spectrum)


def linear_spectrum_02(peptide):
    """
    string -> [int]
    algorithm in textbook
    >>> linear_spectrum_02('NQEL')
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


def linear_spectrum_03(ipeptide):
    """
    [int] -> [int]
    modified version of LinearSpectrum function
    returns linear spectrum from a peptide expressed by a list of integers
    >>> linear_spectrum([113,128,186])
        [0,113,128,186,241,314,427]
    """
    prefix_mass = [0]
    for mass in ipeptide:
        prefix_mass.append(prefix_mass[-1] + mass)
    lspectrum = [0]
    for i in range(len(ipeptide)):
        for j in range(0, len(ipeptide)-i):
            lspectrum.append(prefix_mass[j+i+1] - prefix_mass[j])
    return sorted(lspectrum)


def pretty_print(ls):
    for elm in ls:
        print(elm, end=' ')
    print()


def main():
    try:
        with open('/home/wsl/rosalind/data/ba04j.txt', 'r') as f:
            peptide = f.readline().strip()
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    lspec = linear_spectrum_01(peptide)
    pretty_print(lspec)
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()
