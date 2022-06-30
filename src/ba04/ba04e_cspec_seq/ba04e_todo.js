/*
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
  ║         Peptides <- Expand(Peptides)                         ║ <-- expand
  ║         for each peptide Peptide in Peptides                 ║
  ║             if Mass(Peptide) = ParentMass(Spectrum)          ║
  ║                 if Cyclospectrum(Peptide) = Spectrum         ║
  ║                     output Peptide                           ║
  ║                 remove Peptide from Peptides                 ║
  ║             else if Peptide is not consistent with Spectrum  ║ <-- check consistency
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

References:
- For-each over an array in JavaScript
  https://stackoverflow.com/questions/9329446/for-each-over-an-array-in-javascript
- How to compare arrays in JavaScript?
  https://stackoverflow.com/questions/7837456/how-to-compare-arrays-in-javascript
  JSON.stringify(a1) === JSON.stringify(a2);
  arr1.every((value, index) => value == arr2[index])
- How to find the sum of an array of numbers
  https://stackoverflow.com/questions/1230233/how-to-find-the-sum-of-an-array-of-numbers
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

// helper fx.: array comparison
function arr_comp_loop(arr1, arr2) {
    if (arr1.length !== arr2.length) return false;
    for (var i = 0, len = arr1.length; i < len; i++) {
        if (arr1[i] !== arr2[i]) {
            return false;
        }
    }
    return true;
}

// helper fx.: array comparison
function arr_comp(arr1, arr2) {
    if (arr1.length !== arr2.length) return false;
    return arr1.every((value, index) => value == arr2[index]);
}

// helper fx.: array comparison
function arr_comp(arr1, arr2) {
    return JSON.stringify(arr1) === JSON.stringify(arr2);
}

// linear spectrum of inters, from integer peptides
function lspectrum_from_ipeptide(ipeptide) {
    /**
     * [int] -> [int]
     * return a list of masses
     * >>> lspectrum_from_ipeptide([113,114,128,129])  // LNQE
     *     [0,113,114,128,129,227,242,257,355,371,484]
     */
    let prefixmass = [0];
    for (const mass of ipeptide) {
        prefixmass.push(prefixmass[prefixmass.length - 1] + mass);
    }
    let lspec = [0];
    for (let i = 0; i < ipeptide.length; i++) {
        for (let j = i + 1; j < ipeptide.length + 1; j++) {
            lspec.push(prefixmass[j] - prefixmass[i]);
        }
    }
    return lspec.sort((a, b) => a - b);
}

// cyclo spectrum of inters, from integer peptides
function cspectrum_from_ipeptide(ipeptide) {
    /**
     * [int] -> [int]
     * returns a cyclospectrum from integer peptide
     * >>> cspectrum_from_ipeptide([114,128,129])
     *     [0,114,128,129,242,243,257,371]
     */
    let prefixmass = [0];
    for (const mass of ipeptide) {
        prefixmass.push(prefixmass[prefixmass.length - 1] + mass);
    }
    const mass = prefixmass[prefixmass.length - 1];
    let cspec = [0];
    for (let i = 0; i < ipeptide.length; i++) {
        for (let j = i + 1; j < ipeptide.length + 1; j++) {
            cspec.push(prefixmass[j] - prefixmass[i]);
            if (i > 0 && j < ipeptide.length) {
                cspec.push(mass - (prefixmass[j] - prefixmass[i]));
            }
        }
    }
    return cspec.sort((a, b) => a - b);
}

// expand
function expand(ipeptides) {
    /**
     * [[int]] -> [[int]]
     * returns a expanded integer-peptides
     * simply add each of amino acide mass to each of integer peptide one by one
     * >>> expand([[57]])
     *     [[57,57],[57,71],[57,87],[57,97],[57,99],[57,101],[57,103],[57,113],[57,114],
     *     [57,115],[57,128],[57,129],[57,131],[57,137],[57,147],[57,156],[57,163],[57,186]]
     */
    let result = [];
    for (const pep of ipeptides) {
        for (const mass of unique_masses) {
            result.push([...pep, mass]);
        }
    }
    return result;
}

// consistency check
function is_consistent(ipeptide, cspectrum) {
    /**
     * ([int],[int]) -> bool
     * returns true if linear spectrum of ipeptide is consistent with given cyclospectrum
     * >>> is_consistent([114,128],[0,114,128,129,242,243,257,371])  // LQ, LQE
     *     true
     */
    const lspec = lspectrum_from_ipeptide(ipeptide);
    return lspec.reduce((acc, elem) => cspectrum.includes(elem), true);
}

// BA04E: cyclopeptide sequencing
function cyclopeptide_sequencing(cspectrum) {
    /**
     * [int] -> [[int]]
     * returns a list of integer peptides
     * TODO: too slow (>> 5 minutes, failed to get answer)
     */
    let result = [[]];
    let ipeptides = [[]];
    while (ipeptides.length !== 0) {
        ipeptides = expand(ipeptides);
        for (const ipeptide of ipeptides) {
            const mass = ipeptide.reduce((acc, elem) => acc + elem, 0);
            if (mass === cspectrum[cspectrum.length - 1]) {
                const cspec = cspectrum_from_ipeptide(ipeptide);
                const comp = arr_comp(cspec, cspectrum);
                if (comp) result.push(ipeptide);
                ipeptides = ipeptides.filter(e => !(arr_comp(e, ipeptide)));
            } else if (!(is_consistent(ipeptide, cspectrum))) {
                ipeptides = ipeptides.filter(e => !(arr_comp(e, ipeptide)));
            }
        }
    }
    return result;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04e.txt').toString().split("\n");
    const cspec = lines[0].split(' ').map(Number);

    const startTime = performance.now();
    const ipeptides = cyclopeptide_sequencing(cspec);
    let answer = ipeptides[0];
    for (const ipep of ipeptides.slice(1)) {
        answer += ipep.join('-') + " "
    }
    console.log(answer);
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()