/*
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

// BA04C: cyclospectrum
function spectrum_cyclo(peptide) {
    /**
     * str -> [int]
     * returns a cyclospectrum
     * >>> spectrum_cyclo("NQE")
     */
    let prefixmass = [0];
    for (const acid of peptide) {
        prefixmass.push(prefixmass[prefixmass.length - 1] + weights[acid]);
    }
    const mass = prefixmass[prefixmass.length - 1];
    let cspec = [0];
    for (let i = 0; i < peptide.length; i++) {
        for (let j = i + 1; j < peptide.length + 1; j++) {
            cspec.push(prefixmass[j] - prefixmass[i]);
            if (i > 0 && j < peptide.length) {
                cspec.push(mass - (prefixmass[j] - prefixmass[i]));
            }
        }
    }
    return cspec.sort((a, b) => a - b);
}

// BA04F: score
function score_cyclo_wrong(peptide, cspectrum) {
    /**
     * (str,[int]) -> int
     * returns a score
     * >>> score_cyclo_wrong("NQEL",[0,99,113,114,128,227,257,299,355,356,370,371,484])
     *     11
     */
    const cspec = spectrum_cyclo(peptide);
    return cspec.reduce((acc, e) => (cspectrum.includes(e)) ? acc + 1 : acc, 0);
}

// BA04F: score
function score_cyclo(peptide, cspectrum) {
    /**
     * (str,[int]) -> int
     * returns cyclo score, non-destructive
     */
    let cnt = 0;
    const theoretical_cspectrum = spectrum_cyclo(peptide);
    let [i, j] = [0, 0];
    while (i !== theoretical_cspectrum.length && j !== cspectrum.length) {
        if (theoretical_cspectrum[i] === cspectrum[j]) {
            cnt++; i++; j++;
        } else {
            if (theoretical_cspectrum[i] < cspectrum[j]) i++;
            if (theoretical_cspectrum[i] > cspectrum[j]) j++;
        }
    }
    return cnt;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04f.txt').toString().split("\n");
    const peptide = lines[0];
    const cspectrum = lines[1].split(' ').map(Number);

    const startTime = performance.now();
    console.log(score_cyclo(peptide, cspectrum));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
