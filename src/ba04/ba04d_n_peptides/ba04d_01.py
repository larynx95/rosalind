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

-------------------------------------------------

Dynamic programming:

  (1) initial state of counting dictionary
  0:0, 1:0, ... , 56:0, 57:1, 58:1, ... , 186:1, 187:0, 188:0, ... , 1024:1
  --------(A)---------  ---------(B)-----------  ----------(C)-------------
  always zero           always one               changing part

  (2) accumulating step for counting dictionary
  the index (actually mass) is the sum of masses in part B.
  for example,
    Mass('GG') == 114  --> 114 - 57 = 57 --> Dict[57] == 1 > 0  --> Dict[114] += Dict[114-57]

  (3) get the the number of peptides of given mass

'''

import time


# given info
weight_table ={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,'C':103,'I':113,'L':113,'N':114,
               'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,'R':156,'Y':163,'W':186}
aacids = 'GASPVTCILNDKQEMHFRYW'
aacid_masses  = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
unique_masses = [57,71,87,97,99,101,103,113,    114,115,128,    129,131,137,147,156,163,186]


def CountPeptides(Mass):
    NumPeptides={}
    for i in range(57):
        NumPeptides[i] = 0
    for mass in range(57, Mass + 1):
        NumPeptides[mass] = unique_masses.count(mass)
        for int_mass in unique_masses:
            if mass >= int_mass:
                if NumPeptides[mass - int_mass] > 0:
                    NumPeptides[mass] += NumPeptides[mass - int_mass]
    return NumPeptides[Mass]


def main():
    start_time = time.time()
    print(CountPeptides(1024))
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()


