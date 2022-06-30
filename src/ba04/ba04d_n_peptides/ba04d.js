/*
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
  Does the 'spectrum' mean linear-spectrum or cyclo-spectrum? => Linear

Plan 1.
  * brute-force
  ╔══════════════════════════════════════════════════════════════════════════╗
  ║ Pseudocode:                                                              ║
  ║   BFCYCLOPEPTIDESEQUENCING(Spectrum)                                     ║
  ║       for every Peptide with MASS(Peptide) equal to PARENTMASS(Spectrum) ║
  ║           if Spectrum = CYCLOSPECTRUM(Peptide)                           ║
  ║       output Peptide                                                     ║
  ╚══════════════════════════════════════════════════════════════════════════╝

Plan 2.
  * Dynamic programming

═════════════════════════════════════════════════

References:
-
*/
// #!/usr/bin/env javascript

const weights = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}
const aminos = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
const masses = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186];
const unique_masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186];

// BA04D: compute the number of peptide with given mass
function cout_peptides(mass_given) {
    /**
     * int -> int
     * returns the number of peptides of given mass
     */
    let nums = {};
    // less than 57
    Array(57).fill(0).forEach((e, i) => nums[i] = e);
    // 57 to 186, greater than 186
    for (let m = 57; m < mass_given + 1; m++) {
        nums[m] = unique_masses.includes(m) ? 1 : 0;
        for (const um of unique_masses) {
            if (m >= um && nums[m - um] > 0) {
                nums[m] += nums[m - um];
            }
        }
    }
    return nums[mass_given];
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04d.txt').toString().split("\n");
    const mass_given = lines[0];

    const startTime = performance.now();
    console.log(cout_peptides(mass_given));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()