"""
Rosalind: BA4C
Generate the Theoretical Spectrum of a Cyclic Peptide

The workhorse of peptide sequencing is the mass spectrometer,
an expensive molecular scale that shatters molecules into pieces and then weighs the resulting fragments.
The mass spectrometer measures the mass of a molecule in daltons (Da);
1 Da is approximately equal to the mass of a single nuclear particle (i.e., a proton or neutron).

We will approximate the mass of a molecule by simply adding the number of protons
and neutrons found in the molecule’s constituent atoms,
which yields the molecule’s integer mass.
For example, the amino acid "Gly", which has chemical formula C2H3ON,
has an integer mass of 57, since 2·12 + 3·1 + 1·16 + 1·14 = 57.
Yet 1 Da is not exactly equal to the mass of a proton/neutron,
and we may need to account for different naturally occurring isotopes of each atom when weighing a molecule.
 As a result, amino acids typically have non-integer masses
 (e.g., "Gly" has total mass equal to approximately 57.02 Da);
 for simplicity, however, we will work with the integer mass table given in Figure 1.

  G  A  S  P  V  T   C   I   L   N   D   K   Q   E   M   H   F   R   Y   W
  57 71 87 97 99 101 103 113 113 114 115 128 128 129 131 137 147 156 163 186

The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide),
is the collection of all of the masses of its subpeptides,
in addition to the mass 0 and the mass of the entire peptide.
We will assume that the theoretical spectrum can contain duplicate elements,
as is the case for "NQEL" (shown in Figure 2), where "NQ" and "EL" have the same mass.

    L   N   Q   E   LN  NQ  EL  QE  LNQ ELN QEL NQE NQEL  <-- cyclic
  0 113 114 128 129 227 242 242 257 355 356 370 371 484

Generating Theoretical Spectrum Problem
Generate the theoretical spectrum of a cyclic peptide.

Given: An amino acid string Peptide.

Return: Cyclospectrum(Peptide).

Sample Dataset
LEQN

Sample Output
0 113 114 128 129 227 242 242 257 355 356 370 371 484

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
      -> LinearSpectrum (BA4J)
      -> LinearScore (BA4K)
      -> PREV: the number of Linear peptides (BA4D)
      ↓
    * Cyclo
      -> HERE: Cyclospectrum (BA4C)
      -> NEXT: CycloScore (BA4F)

    [ summary ]

    # Linear vs. Cyclo

    Linear-peptides, spectrum ┌ Peptides -> Spectrum
                              │   ├ Accumulating
                              │   └ Reusing (or Referencing)
                              └ Masses -> Spectrum
    Cyclo-peptides, spectrum  ┌ Peptides -> Spectrum
                              │   ├ Accumulating
                              │   └ Reusing (or Referencing)
                              └ Masses -> Spectrum

    # Chain of fuctions

    linear_subpeptides -> linear_spectrum_from_subpeptides
    cyclo-peptides -> cyclo_spectrum_from_subpeptides


Plan 1.
- peptides, then spectrum

1. Linear-subpeptides, spectrum
1) peptides, then spectrum
- constructing list of subpeptides

  given peptide: NQEL
  i  j  append
  ----------------------- --> generalized rule
  0  0  peptide[0:1] N        outer loop: i in range(len(peptide))
     1  peptide[1:2] Q        inner loop: j in range(len(peptide)-i)
     2  peptide[2:3] E
     3  peptide[3:4] L        subpeptide = peptide[j:i+j+1]
  1  0  peptide[0:2] NQ
     1  peptide[1:3] QE
     2  peptide[2:4] EL
  2  0  peptide[0:3] NQE
     1  peptide[1:4] QEL
  3  0  peptide[0:4] NQEL

- constructing list of subpeptides, masses

    N   Q   E   L   NQ  QE  EL  NQE QEL NQEL
  0 114 128 129 113 242 257 242 371 370 484
  after sorting
  0 113 114 128 129 242 242 257 370 371 484

  Mass(QE) = 128 + 129 = 257
  PrefixMass = 0,114,242,371,484
               values in above list were calculated as belows...
               0,N,  NQ, NQE,NQEL
               0, N        , NQ             , NQE                    , NQEL
               0, 0+Mass(N), Mass(N)+Mass(Q), Mass(N)+Mass(Q)+Mass(E), Mass(N)+Mass(Q)+Mass(E)+Mass(L)
               0, 0+114    , 114+128=242    , 242+129=371            , 371+113=484
               0, 114      , 242            , 371                    , 484

  i  peptides  peptide[i:j]             Mass(peptide)
    --------------------------------------------------------------------------
  0  N Q E L   [0:1] [1:2] [2:3] [3:4]  p[1]-p[0] p[2]-p[1] p[3]-p[2] p[4]-p[3]
  1  NQ QE EL  [0:2] [1:3] [2:4]                  p[2]-p[0] p[3]-p[1] p[4]-p[2]
  2  NQE QEL   [0:3] [1:4]                                  p[3]-p[0] p[4]-p[1]
  3  NQEL      [0:4]                                                  p[4]-p[0]

        subpeptides = []                     spectrum = [0]
        for i=0 to len(p):                   for i=0 for len(p):
            for j=i+1 to len(p)+1:               for j=0 for len(p)-i:
                subpeptides.append(p[0:j])           spectrum.append(Pref[j+i+1]-Pref[j])

2) masses, then spectrum

2. cyclo-subpeptides, spectrum
1) peptides, then spectrum
(1) constructing cyclo-spectrum list of subpeptides - "accumutating" a result list
  focusing on the letter of amino acid
    NQEL  ''                                   common operation (updating result list)
          N Q E L          '' + N,Q,E,L        add the next aacid to each element
          NQ QE EL LN      N+Q Q+E E+L L+N     add the next aacid to each element
          NQE QEL ELN LNQ  NQ+E QE+L EL+N LN+Q add the next aacid to each element
          NQEL             NQE+L               add the next aacid to each element

  pt='NQEL...', len(pt) == n, ls=['',N,Q,E,L,...]  <-- list of peptides with just one amino acid
  i=0 ['',N,Q,E,L]
      ls[-4]+pt[0], ls[-3]+pt[1], ls[-2]+pt[2], ls[-1]+pt[3]
  i=1 ['',N,Q,E,L,NQ,QE,EL,LN]
      ls[-4]+pt[1], ls[-3]+pt[2], ls[-2]+pt[3], ls[-1]+pt[0]
  i=2 ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ]
      ls[-4]+pt[2], ls[-3]+pt[3], ls[-2]+pt[0], ls[-1]+pt[1]
  i=k ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ, ...]
      ls[-n]+pt[(k+0)%n], ls[-n+1]+pt[(k+1)%n], ... , ls[-n+(n-2)]+pt[(k+(n-2))%n] , ls[-n+(n-1)]+pt[(k+(n-1))%n]
  i=n ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ, ...]
      ls+pt

(2) constructing cyclo-spectrum list of subpeptides - "reusing" a string repeatedly ("NQEL")
    L   N   Q   E   LN  NQ  EL  QE  LNQ ELN QEL NQE NQEL  <-- cyclic
  0 113 114 128 129 227 242 242 257 355 356 370 371 484

  i  peptides         peptide[i:j]       generalized rule
  --------------------------------------------------------------------------
  0  N Q E L          [0]                if i+j+1 <= 4, [j:i+j+1]
                      [1]                if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                      [2]
                      [3]
  1  NQ QE EL LN      [0:2]              if i+j+1 <= 4, [j:i+j+1]
                      [1:3]              if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                      [2:4]
                      [3:]    +[:1]
  2  NQE QEL ELN LNQ  [0:3]              if i+j+1 <= 4, [j:i+j+1]
                      [1:4]              if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                      [2:]    +[:1]
                      [3:]    +[:2]
  3  NQEL             [0:4]

2) masses, then spectrum

Plan 2.
- textbook alrogorithms - charging station (textbook pp.211~212)
- creating a "PrefixMass" list, then get spectrum
  PrefixMass('NQEL') = [0,114,242,371,484]
                          N   NQ  NQE NQEL

1. linear
- Pseudocode:
  ╔═══════════════════════════════════════════════════════════════════════╗
  ║ LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)                     ║
  ║     # A. section for constructing PrefixMass list                     ║
  ║     PrefixMass(0) <- 0                                                ║
  ║     for i <- 1 to |Peptide|                                           ║
  ║         for j <- 1 to 20                                              ║
  ║             if AminoAcid(j) =  i-th amino acid in Peptide             ║
  ║                 PrefixMass(i) <- PrefixMass(i - 1) + AminoAcidMass(j) ║
  ║     # B. section for constructing spectrum                            ║
  ║     LinearSpectrum <- a list consisting of the single integer 0       ║
  ║     for i <- 0 to |Peptide| - 1                                       ║
  ║         for j <- i + 1 to |Peptide|                                   ║
  ║             add PrefixMass(j) - PrefixMass(i) to LinearSpectrum       ║
  ║     return sorted list LinearSpectrum                                 ║
  ╚═══════════════════════════════════════════════════════════════════════╝

2. cyclo
- pseudocode (textbook):
  ╔═════════════════════════════════════════════════════════════════════════════════════╗
  ║ CYCLICSPECTRUM(Peptide, AMINOACID, AMINOACIDMASS)                                   ║
  ║     # A. section for constructing PrefixMass list                                   ║
  ║     PREFIXMASS(0) 0                                                                 ║
  ║     for i <- 1 to Peptide                                                           ║
  ║         for j <- 1 to 20                                                            ║
  ║             if AMINOACID(j) = i-th amino acid in Peptide                            ║
  ║                 PREFIXMASS(i) PREFIXMASS(i - 1) + AMINOACIDMASS(j)                  ║
  ║     peptideMass <- PREFIXMASS(|Peptide|)                                            ║
  ║     # B. section for constructing spectrum list                                     ║
  ║     CyclicSpectrum <- a list consisting of the single integer 0                     ║
  ║     for i <- 0 to |Peptide| - 1                                                     ║
  ║         for j <- i + 1 to |Peptide|                                                 ║
  ║             add PREFIXMASS(j) - PREFIXMASS(i) to CyclicSpectrum                     ║
  ║             if i > 0 and j < |Peptide|                                              ║
  ║                 add peptideMass - (PREFIXMASS(j) - PREFIXMASS(i)) to CyclicSpectrum ║
  ║     return sorted list CyclicSpectrum                                               ║
  ╚═════════════════════════════════════════════════════════════════════════════════════╝

═════════════════════════════════════════════════

References:
- Generating all possible combinations of a list, "itertools.combinations" misses some results
  https://stackoverflow.com/questions/17434070/generating-all-possible-combinations-of-a-list-itertools-combinations-misses
- How to generate all permutations of a list?
  https://stackoverflow.com/questions/104420/how-to-generate-all-permutations-of-a-list
- practice complex for-loops
  drawing pyramid ... These practices are not enough!
"""

import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
       'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES = [57, 71, 87, 97, 99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]


def get_mass(peptide):
    """
    string -> int
    returns mass of a peptide
    """
    return sum([WT_TABLE[aacid] for aacid in peptide])


def linear_subpeptides(peptide):
    """
    string -> [string]
    get all linear subpeptides from peptide
    >>> print(linear_subpeptides('NQEL'))
        ['N', 'Q', 'E', 'L', 'NQ', 'QE', 'EL', 'NQE', 'QEL', 'NQEL']
    """
    subpeptides = []
    for i in range(len(peptide)):
        for j in range(len(peptide) - i):
            subpeptides.append(peptide[j:i+j+1])
    return subpeptides


def linear_spectrum_from_subpetides(subpeptides):
    """
    [string] -> [int]
    """
    spectrum = [0]
    for pep in subpeptides:
        spectrum.append(get_mass(pep))
    return sorted(spectrum)


def linear_spectrum(peptide):
    """
    string -> [int]
    returns sorted linear spectrum of a peptide
    >>> linear_spectrum('NQEL')
        [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    """
    subpeptides = []
    for i in range(len(peptide)):
        for j in range(len(peptide) - i):
            subpeptides.append(peptide[j:i+j+1])
    spectrum = [0]
    for pep in subpeptides:
        spectrum.append(get_mass(pep))
    return sorted(spectrum)


def linear_spectrum_textbook(peptide):
    """
    string -> [int]
    algorithm in textbook
    >>> linear_spectrum_textbook('NQEL')
    """
    # constructing cummulative PrefixMass list
    prefix_mass = [0]
    for i in range(len(peptide)):
        for j in range(len(AACIDS)):
            if peptide[i] == AACIDS[j]:
                prefix_mass.append(prefix_mass[-1] + AACID_MASSES[j])
    # constructing spectrum (TODO: This part is more difficult than I thought.)
    linear_spectrum = [0]
    for i in range(len(peptide)):
        for j in range(0, len(peptide)-i):
            linear_spectrum.append(prefix_mass[j+i+1] - prefix_mass[j])
    return sorted(linear_spectrum)


def cyclo_subpeptides(peptide):
    """
    string -> [string]
    get all cyclic subpeptides from peptide
    complex for-loops (not easy, use pencil and paper)
    >>> cyclo_subpeptides('NQEL')
        ['','N','Q','E','L','NQ','QE','EL','LN','NQE','QEL','ELN','LNQ','NQEL']
    """
    ls = ['']
    n = len(peptide)
    for aacid in peptide:
        ls.append(aacid)
    for i in range(len(peptide)-2):
        temp = []
        for j in range(len(peptide)):
            combined = ls[-n+j] + peptide[(i+j+1)%n]
            temp.append(combined)
        ls.extend(temp)
    ls.append(peptide)
    return ls


def cyclo_spectrum_from_subpetides(subpeptides):
    """
    [string] -> [int]
    returns sorted list of peptide weights from list of peptides
    >>> cyclo_spectrum_from_subpetides(cyclo_subpeptides('NQEL'))
        [0,113,114,128,129,227,242,242,257,355,356,370,371,484]
    """
    result = []
    for pep in subpeptides:
        result.append(get_mass(pep))
    return sorted(result)


def cyclo_spectrum(peptide):
    """
    [string] -> [int]
    solution by my Plan 1.
    returns cyclospecturm of a peptide
    >>> cyclo_spectrum('NQEL')
        [0,113,114,128,129,227,242,242,257,355,356,370,371,484]
    """
    ls = ['']
    n = len(peptide)
    for aacid in peptide:
        ls.append(aacid)
    for i in range(len(peptide)-2):
        temp = []
        for j in range(len(peptide)):
            combined = ls[-n+j] + peptide[(i+j+1)%n]
            temp.append(combined)
        ls.extend(temp)
    ls.append(peptide)
    result = []
    for pep in ls:
        result.append(get_mass(pep))
    return sorted(result)


def cyclo_spectrum_textbook(peptide):
    """
    string -> [int]
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


def main():
    f = open('/home/wsl/rosalind/data/ba04c.txt', 'r')
    peptide = f.readline().strip()
    f.close()

    print([0,1,2,3,4][3:9])

    start_time = time.time()
    # print(cyclo_spectrum_textbook(peptide))
    ls =['NQE', 'QNE', 'QEN', 'NEQ', 'ENQ', 'EQN']
    for e in ls:
        print(cyclo_spectrum(e))
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()
