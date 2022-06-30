/*
Rosalind: BA4L
Trim a Peptide Leaderboard

The Trim algorithm, shown below,
sorts all peptides in Leaderboard according to their scores, resulting in a sorted Leaderboard.
Trim> then retains the "top N scoring peptides including ties",
and removes all other peptides from Leaderboard.

╔══════════════════════════════════════════════════════════════════════════════════╗
║ Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)                         ║
║     for j <- 1 to |Leaderboard|                                                  ║
║         Peptide <- j-th peptide in Leaderboard                                   ║
║         LinearScores(j) <- LinearScore(Peptide, Spectrum)                        ║
║     sort Leaderboard according to the decreasing order of scores in LinearScores ║
║     sort LinearScores in decreasing order                                        ║
║     for j <- N + 1 to |Leaderboard|                                              ║
║         if LinearScores(j) < LinearScores(N)                                     ║
║             remove all peptides starting from the j-th peptide from Leaderboard  ║
║         return Leaderboard                                                       ║
║     return Leaderboard                                                           ║
╚══════════════════════════════════════════════════════════════════════════════════╝

Trim Problem
Trim a leaderboard of peptides.

Given: A leaderboard of linear peptides Leaderboard, a linear spectrum Spectrum, and an integer N.

Return: The top N peptides from Leaderboard scored against Spectrum. Remember to use LinearScore.

Sample Dataset
LAST ALST TLLT TQAS
0 71 87 101 113 158 184 188 259 271 372
2

Sample Output
LAST ALST

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
      -> PREV: CycleoPeptideSequencing (BA4E)
      ↓
    * Leaderboard
      -> NEXT: LeaderboardCyclopeptideSequencing (BA4G)
        -> HERE: Trim (BA4L)

Plan 1.
- What's the meaning of the phrase, "top N scoring peptides including ties"?

  sample list: [(10,'a'),(8,'b'),(8,'z'),(3,'c'),(0,'d')]
  top 3 including ties    : [(10,'a'),(8,'b'),(8,'z'),(3,'c')]

  top 1: (10,'a')
  top 2: (8,'b'), (8,'z')
  top 3: (3,'c')

- example
  linear_score('LAST',[0,71,87,101,113,158,184,188,259,271,372])  --> 11
  linear_score('ALST',[0,71,87,101,113,158,184,188,259,271,372])  --> 9
  linear_score('TLLT',[0,71,87,101,113,158,184,188,259,271,372])  --> 5
  linear_score('TQAS',[0,71,87,101,113,158,184,188,259,271,372])  --> 5

  leaderboard = ['LAST','ALST','TLLT','TQAS']
  lspectrum   = [0,71,87,101,113,158,184,188,259,271,372]

  top 1: trim(leaderboard, lspectrum, 1)  --> (11,'LAST')
  top 2: trim(leaderboard, lspectrum, 2)  --> (9,'ALST')
  top 3: trim(leaderboard, lspectrum, 3)  --> (5,'TLLT'), (5,'TQAS')
  top 4: trim(leaderboard, lspectrum, 4)  --> ?

═════════════════════════════════════════════════

References:
- Sorting object property by values
  https://stackoverflow.com/questions/1069666/sorting-object-property-by-values
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
    return lspec.sort((a, b) => a - b);
}

// BA04F: cyclo-score
function score_linear(peptide, lspectrum) {
    /**
     * (str,[int]) -> int
     * returns linear score, non-destructive (BA4K)
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

// BA04L: trim
function trim(lboard, lspectrum, n) {
    /**
     * ([str],[int],int) -> [str]
     * >>> trim(['LAST','ALST','TLLT','TQAS'],[0,71,87,101,113,158,184,188,259,271,372],2)
     *     ['LAST','ALST']
     */
    let tuples = [];
    for (const pep of lboard) {
        tuples.push([pep, score_linear(pep, lspectrum)]);
    }
    tuples.sort((a, b) => b[1] - a[1]);
    let leaders = tuples.map(tup => tup[0]);
    const lscores = tuples.map(tup => tup[1]);
    for (let j = n; j < leaders.length; j++) {
        if (lscores[j] < lscores[n - 1]) {
            leaders = leaders.filter((e, i) => i < j);
            break;
        }
    }
    return leaders;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04l.txt').toString().split("\n");
    const lboard = lines[0].split(/\s/);
    const lspectrum = lines[1].split(/\s/).map(Number);
    const n = parseInt(lines[2]);

    const startTime = performance.now();
    const tops = trim(lboard, lspectrum, n);
    console.log(tops.join(' '));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()