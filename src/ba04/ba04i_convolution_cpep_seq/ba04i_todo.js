/*
Rosalind: BA4I
Implement ConvolutionCyclopeptideSequencing

Given:
An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.

Return:
A cyclic peptide LeaderPeptide with amino acids
taken only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200,
and where the size of Leaderboard is restricted to the top N (and ties).

Sample Dataset
20
60
57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493

Sample Output
99-71-137-57-72-57

═════════════════════════════════════════════════

Info.
  * top N elements with tie?
    [5,5,4,3,3,2,2,1,1]
    top 1 element with tie:  5,5 (tie 5)
    top 2 elements with tie: 5,5
    top 3 elements with tie: 5,5,4
    top 4 elements with tie: 5,5,4,3,3
    top 5 elements with tie: 5,5,4,3,3

  * sample: convolution([0,99,113,114,128,227,257])
    {1:1,14:2,15:2,29:1,30:1,99:2,113:2,114:2,128:2,129:1,143:1,144:1,158:1,227:1,257:1}
    {                        99:2,113:2,114:2,128:2,129:1,143:1,144:1,158:1            }
    top 1 element with tie : 99, 113, 114, 128
    top 2 elements with tie: 99, 113, 114, 128
    top 3 elements with tie: 99, 113, 114, 128
    top 4 elements with tie: 99, 113, 114, 128
    top 5 elements with tie: 99, 113, 114, 128, 129, 143, 144, 158

═════════════════════════════════════════════════

References:
- How To Easily Rank List Without Ties In Excel?
  https://www.extendoffice.com/documents/excel/4796-excel-rank-without-ties.html
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


function expand(ipeptides, masses) {
    /**
     * ([[int]],[int]) -> [[int]]
     * returns a expanded integer-peptides (modifed version)
     * simply add each of amino acide mass to each of integer peptide one by one
     */
    let result = [];
    for (const pep of ipeptides) {
        for (const mass of masses) {
            result.push([...pep, mass]);
        }
    }
    return result;
}


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


function convolution(spectrum) {
    /**
     * [int] -> [int]
     * returns a convolution of spectrum in decreasing order
     */
    spectrum.sort((a, b) => a - b);
    let freq_obj = {};
    for (let i = 1; i < spectrum.length; i++) {
        for (let j = 0; j < i; j++) {
            const diff = spectrum[i] - spectrum[j];
            if (!(diff in freq_obj)) freq_obj[diff] = 1;
            else freq_obj[diff]++;
        }
    }
    let freq_list = [];
    for (const key in freq_obj) {
        freq_list.push([parseInt(key), freq_obj[key]]);
    }
    freq_list.sort((a, b) => b[1] - a[1]);
    let answer = []
    for (const tup of freq_list) {
        if (tup[0] === 0) continue;  // <-- Be careful! No zeros!
        answer.push(...Array(tup[1]).fill(tup[0]));
    }
    return answer;
}


function convolution_masses(spectrum, m) {
    /**
     * ([int],int) -> [int]
     * returns a the list of elements in the convolution of Spectrum in decreasing order of their multiplicities.
     * (only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, with tie)
     * >>> convolution_masses([0,99,113,114,128,227,257],4)  // [99,113,114,128]
     * >>> convolution_masses([0,99,113,114,128,227,257],5)  // [99,113,114,128,129,143,144,158]
     */
    spectrum.sort((a, b) => a - b);
    let freq_obj = {};
    for (let i = 1; i < spectrum.length; i++) {
        for (let j = 0; j < i; j++) {
            const diff = spectrum[i] - spectrum[j];
            if (diff < 57 || diff > 200) continue;  // filtering by mass
            if (!(diff in freq_obj)) freq_obj[diff] = 1;
            else freq_obj[diff]++;
        }
    }
    let freq_list = [];
    for (const key in freq_obj) {
        freq_list.push([parseInt(key), freq_obj[key]]);
    }
    freq_list.sort((a, b) => b[1] - a[1]);
    let masses = [];
    for (const tup of freq_list) {
        if (tup[1] >= freq_list[m - 1][1]) masses.push(tup[0]);
    }
    masses.sort((a, b) => a - b);
    return masses;
}

// BA04I: ConvolutionCyclopeptideSequencing
function cyclopeptide_sequencing_convolution(spectrum, m, n) {
    /**
     * (int,int,[int]) -> [int]
     * returns a integer peptide
     */
    const masses = convolution_masses(spectrum, m);
    let lboard = [[]];
    lpeptide = [];
    while (lboard.length !== 0) {
        lboard = expand(lboard, masses);
        for (const ipep of lboard) {
            const mass = ipep.reduce((a, b) => a + b, 0);
            const parentmass = spectrum[spectrum.length - 1];
            if (mass === parentmass) {
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
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04i.txt').toString().split("\n");
    const m = parseInt(lines[0]);
    const n = parseInt(lines[1]);
    const spectrum = lines[2].split(' ').map(Number);

    const startTime = performance.now();
    const leader_peptide = cyclopeptide_sequencing_convolution(spectrum, m, n);
    console.log(leader_peptide.join('-'));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()