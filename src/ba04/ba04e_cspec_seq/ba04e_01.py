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

═════════════════════════════════════════════════

Preparations.
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

Plan 1. (failed)
- brute-force algorithm
- steps:
  a. ideal experimental spectrum -> a list of single amino acid masses
  b. a list of single amino acid masses -> a list of "ALL" permutations (<-- too many, TODO: Fix this.)
  c. compare theoretical and experimental cyclo-spectrum
  d. final: get a result list

Plan 2.
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

Plan 3.
- using algorithm in textbook
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
      PVC       PVT*      PTP       PTV*      PCV
      99-97-103 99-97-101 99-101-97 99-103-97 101-97-99
      VPC       VPT       VTP*      VCP       TPV
      101-97-103 101-99-97 103-97-101 103-97-99 103-99-97
      TPC        TVP*      CPT        CPV*      CVP


═════════════════════════════════════════════════

References:
- Python | Difference of two lists including duplicates
  https://www.geeksforgeeks.org/python-difference-of-two-lists-including-duplicates/
'''

import time
import math
import sys
import itertools


WT_TABLE ={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
               'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDs = 'GASPVTCILNDKQEMHFRYW'
AACID_MASSES  = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
UNIQUE_MASSES = [57,71,87,97,99,101,103,113,    114,115,128,    129,131,137,147,156,163,186]


def num_subpeptides_linear(len_peptide):
    '''
    returns the number of linear-subpeptides from a n-sized peptide
    >>> num_subpeptides_linear(len('NQEL'))
        11
    '''
    if len_peptide == 1:
        return len_peptide + 1
    else:
        return len_peptide + num_subpeptides_linear(len_peptide - 1)


def num_subpeptides_cyclo(len_peptide):
    '''
    returns the number of cyclo-subpeptides from a n-sized peptide
    >>> num_subpeptides_cyclo(len('NQEL'))
        14
    '''
    if len_peptide == 1:
        return 2
    else:
        return len_peptide * (len_peptide - 1) + 2


def get_answer_polynomial(a,b,c):
    '''
    returns positive integer answer for a polynomial equation 2nd degree
    ax^2 + bx + c = 0
    '''
    answers = ((-b+math.sqrt(b**2-4*a*c))/2*a, (-b-math.sqrt(b**2-4*a*c))/2*a)
    for ans in answers:
        if ans > 0:
            return int(ans)


def linear_spectrum_from_masses(ls):
    '''
    [int] -> [int]
    algorithm in textbook
    >>> linear_spectrum_from_masses([113,128,186])
        [0,113,128,186,241,314,427]
    '''
    prefix_mass = [0]
    for mass in ls:
        prefix_mass.append(prefix_mass[-1] + mass)
    lspectrum = [0]
    for i in range(len(ls)):
        for j in range(0, len(ls)-i):
            lspectrum.append(prefix_mass[j+i+1] - prefix_mass[j])
    return sorted(lspectrum)


def cyclo_spectrum_from_masses(ls):
    '''
    [int] -> [int]
    returns a cyclo-spectrum from a list of single amino acid masses
    >>> cyclo_spectrum_from_masses([113,128,186])
        [0,113,128,186,241,299,314,427]
    '''
    prefix_mass = [0]
    for mass in ls:
        prefix_mass.append(prefix_mass[-1] + mass)
    cspectrum = [0]
    peptide_mass = prefix_mass[-1]
    for i in range(len(ls)):
        for j in range(i+1, len(ls)+1):
            cspectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(ls):
                cspectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cspectrum)


def single_aacid_masses(spectrum):
    '''
    [int] -> [int]
    returns a list of single amino acid masses (both linear and cyclic peptides)
    >>> single_aacid_masses([0,113,128,186,241,299,314,427])
        [113, 128, 186]
    '''
    for i in range(len(spectrum)):
        if spectrum[i] > 186:
            return spectrum[1:i]


def all_perms(ls):
    '''
    [a] -> [[a]]
    >>> list(all_perms('NQE'))
        ['NQE','QNE','QEN','NEQ','ENQ','EQN']
    >>> list(all_perms([113,128,186]))
        [[113,128,186],[128,113,186],[128,186,113],[113,186,128],[186,113,128],[186,128,113]]
    '''
    if len(ls) <=1:
        yield ls
    else:
        for perm in all_perms(ls[1:]):
            for i in range(len(ls)):
                yield perm[:i] + ls[0:1] + perm[i:]


def is_consistent(ls_peptide_masses, spectrum):
    '''
    ([int],[int]) -> Bool
    cspectrum = [0,97,97,99,101,103,196,198,198,200,202,295,297,299,299,301,394,396,398,400,400,497]
    >>> is_consistent([97,101,103], cspectrum)
        False
    >>> is_consistent([97,99,103], cspectrum)
        True
    '''
    lspectrum = linear_spectrum_from_masses(ls_peptide_masses)
    for mass in lspectrum:
        if mass not in spectrum:
            return False
    return True


def is_consistent(ls_peptide_masses, spectrum):
    '''
    cspectrum = [0,97,97,99,101,103,196,198,198,200,202,295,297,299,299,301,394,396,398,400,400,497]
    LinearSpectrum('NQEL') = [0,113,114,128,129,    242,242,257,        370,371,484]
    CycloSpectrum('NQEL') =  [0,113,114,128,129,227,242,242,257,355,356,370,371,484]
    >>> is_consistent([97,101,103], cspectrum)
        False
    >>> is_consistent([97,99,103], cspectrum)
        True
    '''
    lspectrum = linear_spectrum_from_masses(ls_peptide_masses)
    for mass in lspectrum:
        if mass not in spectrum:
            return False
    return True


def is_sum_of_mass_in_spectrum(masses, spectrum):
    '''
    ([int],[int]) -> Bool
    '''
    return sum(masses) in spectrum


def bf_cyclo_peptide_sequencing_ver01_slow(cspectrum):
    '''
    [int] -> [[int]]
    brute-force algorithm: too slow if spectrum is long
    This function only works for short spectrum.
    >>> bf_cyclo_peptide_sequencing_ver01_slow([0,113,128,186,241,299,314,427])
        [[113,128,186],[128,113,186],[128,186,113],[113,186,128],[186,113,128],[186,128,113]]
    '''
    masses = single_aacid_masses(cspectrum)
    perms = all_perms(masses)
    for perm in perms:
        if cyclo_spectrum_from_masses(perm) == cspectrum:
            yield perm


def all_perms_improved(ls, cspectrum):
    if len(ls) <=1:
        yield ls
    else:
        for perm in all_perms_improved(ls[1:], cspectrum):
            for i in range(len(ls)):
                ls = perm[:i] + ls[0:1] + perm[i:]
                if is_consistent(ls, cspectrum):
                    yield ls


def bf_cyclo_peptide_sequencing_ver02_WRONG(cspectrum):
    '''
    [int] -> [[int]]
    returns wrong result, TODO: Fix this function
    >>> bf_cyclo_peptide_sequencing_ver02_WRONG([0,113,128,186,241,299,314,427])
        [[113,128,186],[128,113,186]]
        instead of
        [[113,128,186],[128,113,186],[128,186,113],[113,186,128],[186,113,128],[186,128,113]]
    '''
    masses = single_aacid_masses(cspectrum)
    perms = all_perms_improved(masses, cspectrum)
    for perm in perms:
        if cyclo_spectrum_from_masses(perm) == cspectrum:
            yield perm


def expand_ver01(peptides):
    '''
    [[int]] -> [[int]]
    expand peptide length by one-amino acid (len(result) == len(peptides) * 18)
    >>> expand_ver01([[57]])
        [[57, 57],  [57, 71],  [57, 87],  [57, 97],  [57, 99],  [57, 101],
         [57, 103], [57, 113], [57, 114], [57, 115], [57, 128], [57, 129],
         [57, 131], [57, 137], [57, 147], [57, 156], [57, 163], [57, 186]]
    '''
    new_peptides = []
    for pep in peptides:
        for mass in UNIQUE_MASSES:
            new_peptides.append(pep + [mass])
    return new_peptides


def expand(peptides, spectrum):
    '''
    [[int]] -> [[int]]
    expand peptide length by one-amino acid
    >>> expand([[113]], [0,113,128,186,241,299,314,427])
        [[113, 128], [113, 186]]
    '''
    new_peptides = []
    for pep in peptides:
        for mass in UNIQUE_MASSES:
            temp_expanded = pep + [mass]
            if is_consistent(temp_expanded, spectrum):
                new_peptides.append(temp_expanded)
    return new_peptides


'''
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
'''


def pretty_print(lls):
    for ls in lls:
        ls = map(str, ls)
        print('-'.join(ls), end=' ')
    print()


def main():
    f = open('/home/wsl/rosalind/data/ba04e.txt', 'r')
    cspectrum = [int(elem.strip()) for elem in f.readline().split()]
    f.close()

    sample_cspectrum = [0,113,128,186,241,299,314,427]

    start_time = time.time()
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()


