"""
Rosalind: BA4G
Implement LeaderboardCyclopeptideSequencing

We have thus far worked with theoretical spectra of cyclic peptides,
in which the mass of every subpeptide is given.
This inflexibility presents a practical barrier,
since mass spectrometers generate spectra that are far from ideal
— they are characterized by having both false masses and missing masses.
A false mass is present in the experimental spectrum but absent from the theoretical spectrum;
a missing mass is present in the theoretical spectrum
but absent from the experimental spectrum (see Figure 1).

theoretical:  0    113 114 128 129 227 242 242 257     355 356 370 371 484
experimental: 0 99 113 114 128     227         257 299 355 356 370 371 484

To generalize "Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum"
to handle "noisy" spectra having false and missing masses,
we need to relax the requirement that a candidate peptide’s
theoretical spectrum must match the experimental spectrum exactly,
and instead incorporate a scoring function that will select the peptide
whose theoretical spectrum matches the given experimental spectrum the most closely.
Given a cyclic peptide Peptide and a spectrum Spectrum,
we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum.
Recalling Figure 1, if

Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},
then Score("NQEL", Spectrum) = 11.

To limit the number of candidate peptides under consideration,
we will use a Leaderboard, which holds the N highest scoring candidate peptides for further extension.
At each step, we will expand all candidate peptides found in Leaderboard
by adding every possible amino acid to the end.
Then, we will eliminate those peptides whose newly calculated scores
are not high enough to keep them on the Leaderboard.
This idea is similar to the notion of a "cut" in a golf tournament;
after the cut, only the top N golfers are allowed to play in the next round,
since they are the only players who have a reasonable chance of winning.

To be fair, a cut should include anyone who is tied with the Nth-place competitor.
Thus, Leaderboard should be trimmed down to the "N highest-scoring peptides including ties",
which may include more than N peptides.
Given a list of peptides Leaderboard, a spectrum Spectrum, and an integer N,
Cut(Leaderboard, Spectrum, N)
returns the top N highest-scoring peptides in <Leaderboard (including ties) with respect to Spectrum.

We now introduce LEADERBOARDCYCLOPEPTIDESEQUENCING.
In what follows, the 0-peptide is the peptide "" containing no amino acids.

╔══════════════════════════════════════════════════════════════════════════════╗
║ LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)                               ║
║     Leaderboard <- {0-peptide}                                               ║
║     LeaderPeptide <- 0-peptide                                               ║
║     while Leaderboard is non-empty                                           ║
║         Leaderboard <- Expand(Leaderboard)                                   ║
║         for each Peptide in Leaderboard                                      ║
║             if Mass(Peptide) = ParentMass(Spectrum)                          ║
║                 if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum) ║
║                     LeaderPeptide <- Peptide                                 ║
║             else if Mass(Peptide) > ParentMass(Spectrum)                     ║
║                 remove Peptide from Leaderboard                              ║
║         Leaderboard <- Cut(Leaderboard, Spectrum, N)                         ║
║     output LeaderPeptide                                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

Implement LeaderboardCyclopeptideSequencing
Given: An integer N and a collection of integers Spectrum.

Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).

Sample Dataset
10
0 71 113 129 147 200 218 260 313 331 347 389 460

Sample Output
113-147-71-129

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
      -> CycloScore (BA4F)
      -> CycleoPeptideSequencing (BA4E)
      ↓
    * Leaderboard
      -> HERE: LeaderboardCycloPeptideSequencing (BA4G)
        -> PREV: Trim (BA4L)

Plan 1.
- chains of Functions


  ┌ for-loop (very, very complex version), list slicing, recursion ...
  │
  │ ┌ linear_spectrum ─> linear_score ┌ is_consistent -> expand ┐
  └─┼                                 └ trim ───────────────────┼─┐
    └ cyclo_spectrum ──> cyclo_score ───────────────────────────┘ │
                                                                  │
                          leaderboard_cyclopeptide_sequencing <───┘

- Some functions need to be modified to have an integer list as a factor.

═════════════════════════════════════════════════

References:
- How do I check if a list is empty?
  https://stackoverflow.com/questions/53513/how-do-i-check-if-a-list-is-empty
  [[]] : True
  []   : False
- CAUTION: list can be changed within for-loop
- CAUTION: Python is not strict-typed language.
           It takes unnecessary time to figure out what the return type of the function is.
           Therefore, it is necessary to write down the type of function explicitly.
- CAUTION: Function may not always return the right result.
           In such a case, it is difficult to know which part is wrong.
           TODO: How can I deal with this problem?
"""

#!/usr/bin/env python3
import time


WT_TABLE = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
            'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
AACIDS = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']
AACID_MASSES  = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
UNIQUE_MASSES = [57,71,87,97,99,101,103,113,    114,115,128,    129,131,137,147,156,163,186]


def get_mass(peptide):
    """
    string -> int
    returns mass of a peptide
    """
    return sum([WT_TABLE[aacid] for aacid in peptide])


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


def linear_score(ipeptide, lspectrum):
    """
    ([int],[int]) -> int
    better version, cus lspectrum is intact
    >>> linear_score([114,128,129,113], [0,99,113,114,128,227,257,299,355,356,370,371,484])
        8
    """
    common = []
    theoretical_lspectrum = linear_spectrum(ipeptide)
    for mass in lspectrum:
        if mass in theoretical_lspectrum:
            theoretical_lspectrum.remove(mass)  # <-- theoretical_lspectrum is modifed here, instead of lspectrum
            common.append(mass)
    return len(common)


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


def cyclo_score(ipeptide, cspectrum):
    """
    ([int],[int]) -> int
    >>> cyclo_score([114,128,129,113], [0,99,113,114,128,227,257,299,355,356,370,371,484])
    """
    common = []
    theoretical_cspectrum = cyclo_spectrum(ipeptide)
    for mass in cspectrum:
        if mass in theoretical_cspectrum:
            theoretical_cspectrum.remove(mass)
            common.append(mass)
    return len(common)


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


def trim(ileader_board, cspectrum, n):
    """
    ([int],[int],int) -> set([int])
    modified version of trim - type changed
    >>> trim(['LAST','ALST','TLLT','TQAS'],[0,71,87,101,113,158,184,188,259,271,372],2)
        {'LAST','ALST'}
    >>> trim([[113,71,87,101],[71,113,87,101],[101,113,113,101],[101,128,71,87]],[0,71,87,101,113,158,184,188,259,271,372],2)
        [[113,71,87,101],[71,113,87,101]]
    """
    # get a sorted list of tuples (score, peptide)
    scores = []
    for ipeptide in ileader_board:
        score = linear_score(ipeptide, cspectrum)
        scores.append((score, ipeptide))
    scores.sort(reverse=True)
    # trim
    for j in range(n, len(scores)):
        if scores[j][0] <= scores[n-1][0]:
            scores = scores[:j]
            break
    lboard = []          # can't use set type, unhashable type: 'list'
    for score, ipeptide in scores:
        lboard.append(ipeptide)
    return lboard


def leaderboard_cyclopeptide_sequencing(cspectrum, n):
    """
    ([int],int) -> [int]
    >>> leaderboard_cyclopeptide_sequencing([0,71,113,129,147,200,218,260,313,331,347,389,460],10)
        [129,71,147,113]
    """
    lboard = [[]]
    leader_ipeptide = []
    while lboard:
        lboard = expand(lboard, cspectrum)
        for ipep in lboard:
            if sum(ipep) == cspectrum[-1]:
                if cyclo_score(ipep, cspectrum) > cyclo_score(leader_ipeptide, cspectrum):
                    leader_ipeptide = ipep
            elif sum(ipep) > cspectrum[-1]:
                lboard.remove(ipep)
        lboard = trim(lboard, cspectrum, n)
    return leader_ipeptide


def main():
    f = open('/home/wsl/rosalind/data/ba04g.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    num = int(lines[0])
    cspec = list(map(int, lines[1].split(' ')))
    f.close()

    start_time = time.time()
    result = list(map(str, leaderboard_cyclopeptide_sequencing(cspec, num)))
    print('-'.join(result))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
