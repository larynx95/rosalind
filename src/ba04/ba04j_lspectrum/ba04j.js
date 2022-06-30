/*
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
*/
// #!/usr/bin/env javascript

const weights = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}
const aminos = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
const masses = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

// BA04J: linear spectrum
function spectrum_linear(peptide) {
    /**
     * str -> [int]
     * return a list of masses (BA4J)
     * >>> spectrum_linear("NQEL")
     *     [0,113,114,128,129,242,242,257,370,371,484]
     */
    let prefixmass = [0];
    for (const acid of peptide) {
        prefixmass.push(prefixmass[prefixmass.length - 1] + weights[acid]);
    }
    let lspec = [0];
    for (let i = 0; i < peptide.length; i++) {
        for (let j = i + 1; j < peptide.length + 1; j++) {
            lspec.push(prefixmass[j] - prefixmass[i]);
        }
    }
    return lspec.sort((a, b) => a - b);  // TODO: Be careful!! "sort((a, b) => a - b)", not just "sort()"
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04j.txt').toString().split("\n");
    const peptide = lines[0];

    const startTime = performance.now();
    const lspec = spectrum_linear(peptide);
    console.log(lspec.join(' '));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()