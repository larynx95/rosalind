"""
Rosalind: BA4E
Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

In "Compute the Number of Peptides of Given Total Mass",
we first encountered the problem of reconstructing a cyclic peptide from its theoretical spectrum;
this problem is called the Cyclopeptide Sequencing Problem and is given below.
It is solved by the following algorithm.

  ╔══════════════════════════════════════════════════════════════╗
  ║ CYCLOPEPTIDESEQUENCING(Spectrum)                             ║
  ║     Peptides <- a set containing only the empty peptide      ║
  ║     while Peptides is nonempty                               ║
  ║         Peptides <- Expand(Peptides)                         ║
  ║         for each peptide Peptide in Peptides                 ║
  ║             if Mass(Peptide) = ParentMass(Spectrum)          ║
  ║                 if Cyclospectrum(Peptide) = Spectrum         ║
  ║                     output Peptide                           ║
  ║                 remove Peptide from Peptides                 ║
  ║             else if Peptide is not consistent with Spectrum  ║
  ║                 remove Peptide from Peptides                 ║
  ╚══════════════════════════════════════════════════════════════╝

Cyclopeptide Sequencing Problem
Given an ideal experimental spectrum,
find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

Given: A collection of (possibly repeated) integers Spectrum
corresponding to an ideal experimental spectrum.

Return: Every amino acid string Peptide such that
        Cyclospectrum(Peptide) = Spectrum (if such a string exists).

Sample Dataset
0 113 128 186 241 299 314 427

Sample Output
186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

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
      -> Cyclospectrum (BA4C)
      -> PREV: CycloScore (BA4F)
      -> HERE: CycleoPeptideSequencing (BA4E)
      ↓
    * Leaderboard
      -> LeaderboardCyclopeptideSequencing (BA4G)
        -> NEXT: Trim (BA4L)

Info:
  1. the length of a peptide -> the length of spectrum
    1) linear-spectrum

      peptide  subpeptides (including '', empty peptide)              No.
      ---------------------------------------------------------------------------------
      N        N                                                              1 + 1 = 2
      NQ       N Q        NQ                                                2+1 + 1 = 4
      NQE      N Q E      NQ QE        NQE                                3+2+1 + 1 = 7
      NQEL     N Q E L    NQ QE EL     NQE QEL      NQEL                4+3+2+1 + 1 = 11
      NQELG    N Q E L G  NQ QE EL LG  NQE QEL ELG  NQEL QELG  NQELG  5+4+3+2+1 + 1 = 16
                                                                                num = n! + 1
    2) cyclo-spectrum

      peptide  subpeptides (including '', empty peptide)                                  No.
      ------------------------------------------------------------------------------------------------
      N        N                                                                          2  = 1*0 + 2
      NQ       N Q         NQ                                                             4  = 2*1 + 2 (including zero)
      NQE      N Q E       NQ QE EN        NQE                                            8  = 3*2 + 2
      NQEL     N Q E L     NQ QE EL LN     NQE QEL ELN LNQ      NQEL                      14 = 4*3 + 2
      NQELG    N Q E L G   NQ QE EL LG GN  NQE QEL ELG LGN GNQ  NQEL QELG ELGN LGNQ GNQE  22 = 5*4 + 2
                                                                                        num  = n*(n-1)+2
  2. the length of spectrum -> the length of a peptide
    1) by solving polynomial equation (No. Not this way.)
      linear-spectrum: n! + 1 - len(linear-spectrum) = 0
      cyclo-spectrum: n(n-1) + 2 - len(cyclo-spectrum) = 0
    2) by choosing only those that are less than or equal to the heaviest amino acid mass (186).

  3. branch and bound algorithm
    1) "branch" algorithm
      - branching with all 18 amino acids
      - branching with 17 amino acids
      - branching only with peptides which are "consistent"
    2) bound algorithm
      - "consistent" peptide: all masses in LinearSpectrum(peptide) are in the given cyclospectrum
      - compare total mass
      - compare cyclospectrum

═════════════════════════════════════════════════

Plan 1.
- using algorithm in textbook
- functionality:
  a. obtaining cyclospectrum from peptide expressed in a list of integers
  b. obtaining linearspectrum from peptide expressed in a list of integers
  c. expanding peptide length
    - len(int_peptide) * 18
    - len(int_peptide) * 17
    - expand each peptide which is consistent with cyclospectrum
  d. comparing peptide mass with cyclospectrum[-1]
    - for every peptide
    - for peptides with the same number as single amino acids in cyclospectrum (<= 186)
  3. comparing two cyclospectrums
- two parts:
  a. "branch" (expanding)
  b. "bound"
    - LinearSpectrum(a list of peptide masses) vs. cyclo-spectrum
    - ParentMass (sum of all peptide masses in a list) vs. cyclo-spectrum[-1]
- example in textbook
  cyclo-spectrum: 0 97 97 99 101 103 196 198 198 200 202 295 297 299 299 301 394 396 398 400 400 497

  a. initial single amino acids

    97  99 101 103    --> expand --> 18*4=72 peptides --> trimming process (how??)
     P   V   T   C                   (or --> 17*4=68 peptides)

  b. next expanding and trimming

    - itertools.permutation('PVTC', 2):
      PV     PT     PC     VP     VT     VC     TP     TV     CP     CV     trimmed: CT,TC

    - process:
      LinearSpectrum([97,99])  : [0,97,99,196]         PV, selected
      LinearSpectrum([101,103]): [0,101,103,204]  -->  TC, 204 is not in cyclo-spectrum, so this is trimmed
      LinearSpectrum([103,101]): [0,101,103,204]  -->  CT, 204 is not in cyclo-spectrum, so this is trimmed

    - result:
      97-99  97-101 97-103 99-97  99-101 99-103 101-97 101-99 103-97 103-99
      PV     PT     PC     VP     VT     VC     TP     TV     CP     CV

  c. next expanding and trimming

    - itertools.permutations('PVCT', 3):
       selected: PVC PVT PTV PCV VPC VPT VTP VCP TPV TPC TVP CPT CPV CVP
       trimmed : PTC PCT VTC VCT TVC TCP TCV CVT CTP CTV
       new     : PTP

    - process:
      LinearPectrum([97,99,103]) : [0,97,99,103,196,202,299]  --> PVC, selected
      LinearPectrum([97,101,103]): [0,97,101,103,198,204,301] --> PTC, 204 is not in cyclo-spectrum, so this is trimmed

    - result:
      97-99-103 97-99-101 97-101-97 97-101-99 97-103-99     (*) will be trimmed in next iteration (why?)
      PVC       PVT*      PTP       PTV*      PCV           cuz its linearspectrum is not consistent
      99-97-103 99-97-101 99-101-97 99-103-97 101-97-99
      VPC       VPT       VTP*      VCP       TPV
      101-97-103 101-99-97 103-97-101 103-97-99 103-99-97
      TPC        TVP*      CPT        CPV*      CVP

Plan 2. (failed)
- brute-force algorithm
- steps:
  a. ideal experimental spectrum -> a list of single amino acid masses
  b. a list of single amino acid masses -> a list of "ALL" permutations (<-- too many, TODO: Fix this.)
  c. compare theoretical and experimental cyclo-spectrum
  d. final: get a result list

Plan 3.
- Notice! This is wrong approach.
- still brute-force algorithm, but using filtering
- steps:
  a. ideal experimental spectrum -> a list of single amino acid masses
  b. a list of single amino acid masses -> a list of "SELECTED" permutations
    select permutations satisfying LinearSpectrum(perm) is "CONSISTENT" with the experimental spectrum
    discard useless candidates as much as possible
    but ... The attempt to solve the problem after finding all possible permutations may be wrong.
    Let's take a look at the example below.

    97  99 101 103
     P   V   T   C

    97-99  97-101 97-103 99-97  99-101 99-103 101-97 101-99 103-97 103-99
      PV     PT     PC     VP     VT     VC     TP     TV     CP     CV
    These are not permutations at all!

    So attempts to solve problems by permutation are wrong approaches!
    Unfortunately, a lot of time was wasted.
    TODO: Can I expand peptides by permutation way?
          How can I implement function for proper permutations?
  c. compare theoretical and experimental cyclo-spectrum
  d. final: get a result list

═════════════════════════════════════════════════

References:
- mutation of list in each iteration problem
  How to deep copy a list?
  https://stackoverflow.com/questions/17873384/how-to-deep-copy-a-list
  Python List copy()
  https://www.programiz.com/python-programming/methods/list/copy
"""

import time


WT_DICT ={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
            'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = 'GASPVTCILNDKQEMHFRYW'
AACID_MASSES  = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
UNIQUE_MASSES = [57,71,87,97,99,101,103,113,    114,115,128,    129,131,137,147,156,163,186]


def linear_spectrum(ipeptide):
    """
    [int] -> [int]
    algorithm in textbook
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


def cyclo_spectrum(ipeptide):
    """
    [int] -> [int]
    returns a cyclo-spectrum from a list of single amino acid masses
    >>> cyclo_spectrum([113,128,186])
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


def is_consistent(ipeptide, cspectrum):
    """
    ([int],[int]) -> Bool
    cspectrum = [0,97,97,99,101,103,196,198,198,200,202,295,297,299,299,301,394,396,398,400,400,497]
    >>> is_consistent([97,101,103], cspectrum)
        False
    >>> is_consistent([97,99,103], cspectrum)
        True
    """
    for mass in linear_spectrum(ipeptide):
        # TODO: Why should I use linear_spectrum rather than cyclo_spectrum fx here?
        #       This is very important question!
        #       What if I use cyclo_spectrum fx here?
        #       반드시 생각해볼 것!
        if not mass in cspectrum:
            return False
    return True


def expand(ipeptides, cspectrum):
    """
    ([[int]],[int]) -> [[int]]
    >>> cspect = cyclo_spectrum([113,128,186])  # [0, 113, 128, 186, 241, 299, 314, 427]
    >>> expand([[113]], cspect)
        [[113, 128], [113, 186]]
    """
    for peptide in ipeptides.copy():
        ipeptides.remove(peptide)
        for mass in UNIQUE_MASSES:
            new_peptide =  peptide + [mass]
            if is_consistent(new_peptide, cspectrum):
                ipeptides.append(new_peptide)
    return ipeptides


def cyclopeptide_sequencing(cspectrum):
    """
    [int] -> [[int]]
    >>> sample_cspec = [0,113,128,186,241,299,314,427]
    >>> cyclopeptide_sequencing(sample_cspec)
        [[],[113,128,186],[113,186,128],[128,113,186],[128,186,113],[186,113,128],[186,128,113]]
    """
    result = [[]]
    ipeptides = [[]]
    parent_mass = cspectrum[-1]
    while ipeptides:
        ipeptides = expand(ipeptides, cspectrum)
        for ipeptide in ipeptides.copy():
            if sum(ipeptide) == parent_mass:
                if cyclo_spectrum(ipeptide) == cspectrum:
                    result.append(ipeptide)
                ipeptides.remove(ipeptide)
            elif not is_consistent(ipeptide, cspectrum):
                ipeptides.remove(ipeptide)
    return result


def pretty_print(lls):
    for ls in lls:
        ls = map(str, ls)
        print('-'.join(ls), end=' ')
    print()


def main():
    f = open('/home/wsl/rosalind/data/ba04e.txt', 'r')
    cspectrum = [int(elem.strip()) for elem in f.readline().split()]
    f.close()

    sample_cspec = [0,113,128,186,241,299,314,427]

    start_time = time.time()
    #result = cyclopeptide_sequencing(sample_cspec)
    result = cyclopeptide_sequencing(cspectrum)
    pretty_print(result)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
