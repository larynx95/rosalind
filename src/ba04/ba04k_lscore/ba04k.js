/*
Rosalind: BA4K
Compute the Score of a Linear Peptide

Linear Peptide Scoring Problem
Compute the score of a linear peptide with respect to a spectrum.

Given: An amino acid string Peptide and a collection of integers LinearSpectrum.

Return: The linear score of Peptide against Spectrum, LinearScore(Peptide, Spectrum).

Sample Dataset
NQEL
0 99 113 114 128 227 257 299 355 356 370 371 484

Sample Output
8

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
      -> PREV: LinearSpectrum (BA4J)
      -> HERE: LinearScore (BA4K)
      -> NEXT: the number of Linear peptides (BA4D)

Info.
  * two lists
    theoretical : [0    113 114 128 129     242 242 257             370 371 484]
    experimental: [0 99 113 114 128     227         257 299 355 356 370 371 484]
                   ^    ^   ^   ^                   ^               ^   ^   ^
  * problems:
    a. There're duplicates in one or both lists!
       => element should be removed during comparison
       0112 3445 -> *112 3445 -> *11* 3445 -> *11* 3*45 -> *11* 3*4*
       0  22 4 5    *  22 4 5    *  *2 4 5    *  *2 * 5    *  *2 * *
       ^               ^               ^              ^

    b. One of the lists can be modified during iteration!
       => select theoretical linear-spectrum as modifying one
       => TODO: I found weird thing. (below) Why???

═════════════════════════════════════════════════

References:
- Simplest code for array intersection in javascript
  https://stackoverflow.com/questions/1885557/simplest-code-for-array-intersection-in-javascript
- How to sort an array of integers correctly (TODO: Don't miss this! It may save hours even days.)
  https://stackoverflow.com/questions/1063007/how-to-sort-an-array-of-integers-correctly
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
    return lspec.sort((a, b) => a - b);  // "sort()" sorts only strings. Be careful.
}

//
function score_linear_reduce_wrong1(peptide, lspectrum) {
    /**
     * (str,[int]) -> int
     * returns a score
     */
    const lspec = spectrum_linear(peptide);
    return lspec.reduce((acc, e) => (lspectrum.includes(e)) ? acc + 1 : acc, 0);  // reduce() fx wrong! TODO: Why?
}

//
function score_linear_reduce_wrong2(peptide, lspectrum) {
    /**
     * (str,[int]) -> int
     * returns a score
     */
    const lspec = spectrum_linear(peptide);
    return lspec.reduce((a, b) => a + lspectrum.includes(b), 0);  // reduce() fx wrong! TODO: Why?
}

//
function score_linear_filter_wrong(peptide, lspectrum) {
    const lspec = spectrum_linear(peptide);
    return lspec.filter(e => lspectrum.includes(e)).length;
}

// BA04K: linear score
function score_linear(peptide, lspectrum) {
    /**
     * (str,[int]) -> int
     * returns linear score, non-destructive
     * >>> var lspec = spectrum_linear("LIQE");  // [0,113,113,128,129,226,241,257,354,370,483]
     * >>> score_linear(lspec, [0,113,128,128,226,226,241,242,257])
     *     6
     */
    let cnt = 0;
    const theoretical_lspectrum = spectrum_linear(peptide);
    let [i, j] = [0, 0];
    while (i !== theoretical_lspectrum.length && j !== lspectrum.length) {
        if (theoretical_lspectrum[i] === lspectrum[j]) {
            cnt++; i++; j++;
        } else {
            if (theoretical_lspectrum[i] < lspectrum[j]) i++;
            if (theoretical_lspectrum[i] > lspectrum[j]) j++;
        }
    }
    return cnt;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04k.txt').toString().split("\n");
    const peptide = lines[0];
    const lspectrum = lines[1].split(' ').map(Number);

    const startTime = performance.now();
    console.log(score_linear(peptide, lspectrum));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()
