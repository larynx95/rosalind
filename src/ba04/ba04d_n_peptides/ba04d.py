"""
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

Denote the total mass of an amino acid string Peptide as "Mass(Peptide)".
In mass spectrometry experiments, whereas the peptide that generated a spectrum is unknown,
the peptide's mass is typically known and is denoted "ParentMass(Spectrum)".
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
      -> PREV: LinearScore (BA4K)
      -> HERE: the number of Linear peptides (BA4D)
      ↓
    * Cyclo
      -> NEXT: Cyclospectrum (BA4C)

Info:
- weight table of amino acids
  G  A  S  P  V  T   C   I   L   N   D   K   Q   E   M   H   F   R   Y   W
  57 71 87 97 99 101 103 113 113 114 115 128 128 129 131 137 147 156 163 186

- ambiguity of this problem:
  (1) Does the 'spectrum' mean linear-spectrum or cyclo-spectrum? => Linear


Plan 1.
- brute force algorithm
- chain of functions

  parent_mass      ┐
  all_combs_gen    ┤
  all_perms        ┼──> bf_linear_peptides_searching (& counting)
  linear_spectrum  ┘    (It sucks! Don't try this approach.)

  parent_mass      ┐
  all_combs_gen    ┤
  all_perms        ┼──> bf_cyclo_peptides_searching (& counting)
  cyclo_spectrum   ┘

- steps:
  (1) get all combinations of amino acids in given peptide
  (2) get all permutations (peptides) for each combination above
  (3) get linear spectrum of each peptide in above step
  (4) comparison - theoretical vs. experimental
  (5) select peptides with the same linear spectrum as fx argument (experimental)
  spectrum       <-- linear spectrum

- First step: combinations
  Combinations: get all combination of amino acids from the given peptides
  nCr = n! / r! (n - r)!
  total No. = 3C1 + 3C2 + 3C3 = 3! / 2! + 3! / 2! + 3! / 3! = 7
  "NQE"  '' N Q E NQ QE NE NQE
  Problem can be solve only with combination!
  Notice that the result spectrums are not the same!
  The order of connecting amino acids must be considered.
  So permutation required afterward.
    linear_spectrum('NQE')  --> [0,114,128,129,242,257,371]
    linear_spectrum('EQN')  --> [0,114,128,129,242,257,371]
    linear_spectrum('QNE')  --> [0,114,128,129,242,257,371]
    linear_spectrum('NEQ')  --> [0,114,128,129,243,257,371]
    linear_spectrum('ENQ')  --> [0,114,128,129,242,243,371]

- Second step: permutations
  Permutations (without repetition):
  nPr = n! / (n - r)!
  3P3 = 3! / (3 - 3)! = 3! = 6
  combs  perms  linear spectrum
  -----------------------------------------------
  N      N      [0, 114]
  Q      Q      [0, 128]
  E      E      [0, 129]              lots of peptides are meaningless!!
  NQ     NQ     [0, 114, 128, 242]
         QN     [0, 114, 128, 242]
  QE     QE     [0, 128, 129, 257]
         EQ     [0, 128, 129, 257]
  NE     NE     [0, 114, 129, 243]
         EN     [0, 114, 129, 243]
  NQE    NQE    [0, 114, 128, 129, 242, 257, 371]
         NEQ    [0, 114, 128, 129, 243, 257, 371]
         QNE    [0, 114, 128, 129, 242, 243, 371]
         QEN    [0, 114, 128, 129, 243, 257, 371]
         ENQ    [0, 114, 128, 129, 242, 243, 371]
         EQN    [0, 114, 128, 129, 242, 257, 371]

Plan 2.
- discarding meaningless peptides
  ╔══════════════════════════════════════════════════════════════════════════╗
  ║ Pseudocode:                                                              ║
  ║   BFCYCLOPEPTIDESEQUENCING(Spectrum)                                     ║
  ║       for every Peptide with MASS(Peptide) equal to PARENTMASS(Spectrum) ║
  ║           if Spectrum = CYCLOSPECTRUM(Peptide)                           ║
  ║       output Peptide                                                     ║
  ╚══════════════════════════════════════════════════════════════════════════╝

Plan 3.
- Dynamic programming:
  TODO: difficult! Pratice this algorithm!

  (1) initial state of counting dictionary
  0:0, 1:0, ... , 56:0, 57:1, 58:1, ... , 186:1, 187:0, 188:0, ... , 1024:1
  --------(A)---------  ---------(B)-----------  ----------(C)-------------
  always zero           always one               changing part

  (2) accumulating step for counting dictionary
  the index (actually mass) is the sum of masses in part B.
  for example,
    Mass('GG') == 114  --> 114 - 57 = 57 --> Dict[57] == 1 > 0  --> Dict[114] += Dict[114-57]

  (3) get the the number of peptides of given mass

═════════════════════════════════════════════════

References:
- How to generate all permutations of a list?
  https://stackoverflow.com/questions/104420/how-to-generate-all-permutations-of-a-list
- How to get all possible combinations of a list's elements?
  https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
- Where can I find source code for itertools.combinations() function
  https://stackoverflow.com/questions/5731505/where-can-i-find-source-code-for-itertools-combinations-function
"""

from random import sample
import time
import itertools


# given info
weight_table ={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
               'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
aacids = 'GASPVTCILNDKQEMHFRYW'
aacid_masses  = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
unique_masses = [57,71,87,97,99,101,103,113,    114,115,128,    129,131,137,147,156,163,186]


def get_mass(peptide):
    """
    implementation of 'Mass(peptide)' function
    returns mass of a peptide
    """
    return sum([weight_table[aacid] for aacid in peptide])


def parent_mass(spectrum):
    """
    implementation of 'ParentMass(spectrum)' function
    returns the maximum value in the given spectrum
    (element at the last index)
    """
    return spectrum[-1]


def all_combs_rec(text):
    """
    >>> all_combs_rec('NQE')
        ['', 'N', 'Q', 'QN', 'E', 'EN', 'EQ', 'EQN']
    """
    if len(text) == 0:
        return ['']
    cs = []
    for c in all_combs_rec(text[1:]):
        cs += [c, c+text[0]]
    return cs


def all_combs_gen(text):
    """
    all combinations, generator version
    >>> list(all_combs('NQE'))
        ['', 'N', 'Q', 'QN', 'E', 'EN', 'EQ', 'EQN']
    >>> print(*all_combs_gen('LQNE'))
        L Q QL N NL NQ NQL E EL EQ EQL EN ENL ENQ ENQL
    >>> print(*all_combs_gen('NQEL'))
        N Q QN E EN EQ EQN L LN LQ LQN LE LEN LEQ LEQN
    """
    if len(text) == 0:
        yield ''
    else:
        for c in all_combs_gen(text[1:]):
            yield c
            yield c+text[0]


def all_perms(text):
    """
    >>> list(all_perms('NQE'))
        ['NQE', 'QNE', 'QEN', 'NEQ', 'ENQ', 'EQN']
    """
    if len(text) <=1:
        yield text
    else:
        for perm in all_perms(text[1:]):
            for i in range(len(text)):
                yield perm[:i] + text[0:1] + perm[i:]


def linear_spectrum(peptide):
    """
    returns sorted linear spectrum of a peptide
    >>> linear_spectrum('NQEL')
        [0,113,114,128,129,242,242,257,370,371,484]
    >>> linear_spectrum('LQNE')
        [0,113,114,128,129,241,242,243,355,371,484]
    """
    subpeptides = []
    for i in range(len(peptide)):
        for j in range(len(peptide) - i):
            subpeptides.append(peptide[j:i+j+1])
    spectrum = [0]
    for pep in subpeptides:
        spectrum.append(get_mass(pep))
    return sorted(spectrum)


def cyclo_spectrum(peptide):
    """
    algorithm in textbook
    """
    prefix_mass = [0]
    for aacid in peptide:
        prefix_mass.append(prefix_mass[-1] + weight_table[aacid])
    cspectrum = [0]
    peptide_mass = prefix_mass[-1]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            cspectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cspectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cspectrum)


def bf_linear_peptides_searching(spectrum):
    """
    returns a list of peptides with given linear spectrum
    >>> sample_lspectrum = linear_spectrum('NQE')
    >>> bf_linear_peptides_searching(sample_lspectrum)
        ['EKN', 'NKE', 'EQN', 'NQE']
    """
    found = []
    combs = all_combs_gen(aacids)    # with full set of amino acids... maybe slow
    # combs = all_combs_gen('NQEL')  # <-- this line for testing
    peps = []
    for com in combs:
        if get_mass(com) == spectrum[-1]:  # discard meaningless peptides
            ls_peps = all_perms(com)
            for pep in ls_peps:
                peps.append(pep)
    for pep in peps:
        lspectrum = linear_spectrum(pep)
        if lspectrum == spectrum:
            found.append(pep)
    return found


def bf_cyclo_peptides_searching(spectrum):
    """
    returns a list of peptides with given cyclo spectrum
    >>> sample_lspectrum = cyclo_spectrum('NQE')
    >>> bf_cyclo_peptides_searching(sample_lspectrum)
        ['EKN', 'NKE', 'EQN', 'NQE']
    """
    found = []
    combs = all_combs_gen(aacids)    # with full set of amino acids... maybe slow
    # combs = all_combs_gen('NQEL')  # <-- this line for testing
    peps = []
    for com in combs:
        if get_mass(com) == spectrum[-1]:  # discard meaningless peptides
            ls_peps = all_perms(com)
            for pep in ls_peps:
                peps.append(pep)
    for pep in peps:
        lspectrum = cyclo_spectrum(pep)
        if lspectrum == spectrum:
            found.append(pep)
    return found


def bf_linear_peptides_counting(spectrum):
    """
    >>> bf_linear_peptides_counting([0, 114, 128, 129, 242, 257, 371])
        ['EKN', 'NKE', 'EQN', 'NQE']
    """
    combs = all_combs_gen(aacids)
    peps = []
    for com in combs:
        if get_mass(com) == spectrum[-1]:
            ls_peps = all_perms(com)
            for pep in ls_peps:
                peps.append(pep)
    count = 0
    for pep in peps:
        lspectrum = linear_spectrum(pep)
        if lspectrum == spectrum:
            count += 1
    return count


def bf_counting_peptides_with_given_mass(mass):
    """
    DO NOT execute this function (too slow)
    """
    combs = all_combs_gen(aacids)
    peps = []
    for com in combs:
        if get_mass(com) == mass:
            ls_peps = all_perms(com)
            for pep in ls_peps:
                peps.append(pep)
    count = 0
    for pep in peps:
        if get_mass(pep) == mass:
            count += 1
    return count


def count_peptides(mass_given):
    """
    dynamic programming
    """
    num_peptides = {}
    for i in range(57):
        num_peptides[i] = 0
    for mass in range(57, mass_given + 1):
        num_peptides[mass] = unique_masses.count(mass)
        for int_mass in unique_masses:
            if mass >= int_mass:
                if num_peptides[mass - int_mass] > 0:
                    num_peptides[mass] += num_peptides[mass - int_mass]
    return num_peptides[mass_given]


def main():
    f = open('/home/wsl/rosalind/data/ba04d.txt', 'r')
    mass = int(f.readline().strip())
    f.close()

    # sample specturm for simple testing
    #sample_lspectrum = linear_spectrum('NQE')
    #print(sample_lspectrum)

    # time check 1
    #start_time = time.time()
    #peptides = brute_force_peptide_searching(sample_lspectrum)
    #print(peptides)                # ['EKN', 'NKE', 'EQN', 'NQE']
    # print(linear_spectrum('EKN'))  # [0, 114, 128, 129, 242, 257, 371]
    # print(linear_spectrum('NKE'))  # [0, 114, 128, 129, 242, 257, 371]
    # print(linear_spectrum('EQN'))  # [0, 114, 128, 129, 242, 257, 371]
    # print(linear_spectrum('NQE'))  # [0, 114, 128, 129, 242, 257, 371]
    #print("--- %s seconds ---" % (time.time() - start_time))  # --- 1.1843180656433105 seconds ---

    #start_time = time.time()
    #print(brute_force_peptide_counting(sample_lspectrum))
    #print("--- %s seconds ---" % (time.time() - start_time))  # --- 1.0633888244628906 seconds ---

    #start_time = time.time()
    #print(bf_counting_peptides_with_given_mass(371)) # 336
    #print("--- %s seconds ---" % (time.time() - start_time))  # --- 1.0731542110443115 seconds ---

    start_time = time.time()
    print(count_peptides(1024))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
