/*
Rosalind: BA4G
Implement LeaderboardCyclopeptideSequencing

We have thus far worked with theoretical spectra of cyclic peptides,
in which the mass of every subpeptide is given.
This inflexibility presents a practical barrier,
since mass spectrometers generate spectra that are far from ideal
— they are characterized by having both "false" masses and "missing" masses.
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
we define "Score(Peptide, Spectrum) "as the number of masses shared between "Cyclospectrum(Peptide)" and "Spectrum".
Recalling Figure 1, if

Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},
then
Score("NQEL", Spectrum) = 11.

To limit the number of candidate peptides under consideration,
we will use a Leaderboard, which holds the N highest scoring candidate peptides for further extension.
At each step, we will expand all candidate peptides found in Leaderboard
by adding every possible amino acid to the end.
Then, we will eliminate those peptides whose newly calculated scores
are not high enough to keep them on the Leaderboard.
This idea is similar to the notion of a "cut" in a golf tournament;
after the cut, only the top N golfers are allowed to play in the next round,
since they are the only players who have a reasonable chance of winning.

To be fair, a "cut" should include anyone who is tied with the Nth-place competitor.
Thus, Leaderboard should be trimmed down to the "N highest-scoring peptides including ties",
which may include more than N peptides.
Given a list of peptides Leaderboard, a spectrum Spectrum, and an integer N,
"Cut(Leaderboard, Spectrum, N)"
returns the top N highest-scoring peptides in <Leaderboard (including ties) with respect to Spectrum.

We now introduce LEADERBOARDCYCLOPEPTIDESEQUENCING.
In what follows, the 0-peptide is the peptide "" containing no amino acids.

╔══════════════════════════════════════════════════════════════════════════════╗
║ LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)                               ║
║     Leaderboard <- set containing only the empty peptide                     ║
║     LeaderPeptide <- 0-peptide                                               ║
║     while Leaderboard is non-empty                                           ║
║         Leaderboard <- Expand(Leaderboard)                                   ║ <-- expand
║         for each Peptide in Leaderboard                                      ║
║             if Mass(Peptide) = ParentMass(Spectrum)                          ║ <-- sum of array
║                 if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum) ║ <-- score
║                     LeaderPeptide <- Peptide                                 ║
║             else if Mass(Peptide) > ParentMass(Spectrum)                     ║
║                 remove Peptide from Leaderboard                              ║
║         Leaderboard <- Cut(Leaderboard, Spectrum, N)                         ║ <-- Trim
║     output LeaderPeptide                                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

Implement LeaderboardCyclopeptideSequencing
Given: An integer N and a collection of integers Spectrum.

Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).

Sample Dataset
10
0 71 113 129 147 200 218 260 313 331 347 389 460   // A-I/L-E-F

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

// BA04J: linearspectrum from peptides (integers)
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

// BA04F: cycloscore
function score_linear(ipeptide, lspectrum) {
    /**
     * (str,[int]) -> int
     * returns linear score, non-destructive (BA4K)
     * >>> var lspec = spectrum_linear("LIQE");  // [0,113,113,128,129,226,241,257,354,370,483]
     * >>> score_linear(lspec, [0,113,128,128,226,226,241,242,257])
     *     6
     */
    let cnt = 0;
    const theoretical_lspectrum = lspectrum_from_ipeptide(ipeptide);
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

// BA04C: cyclospectrum, from peptides (intergers)
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

// BA04F: cycleo-score
function score_cyclo(ipeptide, cspectrum) {
    /**
     * (str,[int]) -> int
     * returns cyclo score, non-destructive
     */
    let cnt = 0;
    const theoretical_cspectrum = cspectrum_from_ipeptide(ipeptide);
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

// trim
function trim(lboard, spectrum, n) {
    /**
     * ([str],[int],int) -> [str]
     * >>> trim(['LAST','ALST','TLLT','TQAS'],[0,71,87,101,113,158,184,188,259,271,372],2)
     *     ['LAST','ALST']
     */
    let tuples = [];
    for (const pep of lboard) {
        tuples.push([pep, score_linear(pep, spectrum)]);
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

// expand peptide
function expand(ipeptides, masses) {
    /**
     * [[int]] -> [[int]]
     * returns a expanded integer-peptides
     * simply add each of amino acide mass to each of integer peptide one by one
     * >>> expand([[57]], unique_masses)
     *     [[57,57],[57,71],[57,87],[57,97],[57,99],[57,101],[57,103],[57,113],[57,114],
     *     [57,115],[57,128],[57,129],[57,131],[57,137],[57,147],[57,156],[57,163],[57,186]]
     */
    let result = [];
    for (const pep of ipeptides) {
        for (const mass of masses) {
            result.push([...pep, mass]);
        }
    }
    return result;
}

// BA04G: leaderboard cycleopeptide sequencing
function lb_seq_cyclopeptids(spectrum, n) {
    /**
     * ([int],int) -> str
     * returns a leader peptide
     */
    let lboard = [[]];
    let lpeptide = [];
    while (lboard.length !== 0) {
        lboard = expand(lboard, unique_masses);
        for (const ipep of lboard) {
            const mass = ipep.reduce((a, b) => a + b, 0);
            const parentmass = spectrum[spectrum.length - 1];
            if (mass === parentmass) {  // <-- more verbose than Python
                if (score_cyclo(ipep, spectrum) > score_cyclo(lpeptide, spectrum)) {
                    lpeptide = ipep;
                }
            } else if (mass > parentmass) {
                lboard = lboard.filter(e => e !== ipep);
            }
        }
        lboard = trim(lboard, spectrum, n);
    }
    return lpeptide;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04g.txt').toString().split("\n");
    const n = parseInt(lines[0]);
    const spectrum = lines[1].split(' ').map(Number);

    const startTime = performance.now();
    const leader_peptide = lb_seq_cyclopeptids(spectrum, n);
    console.log(leader_peptide.join('-'));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()