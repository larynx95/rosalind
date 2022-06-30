/*
Rosalind: BA4H
Generate the Convolution of a Spectrum

We define the convolution of a cyclic spectrum by taking all positive differences of masses in the spectrum.
Figure 1 shows the convolution of the theoretical spectrum of "NQEL".

        |  "" L   N   Q   E   LN  NQ  EL  QE  LNQ ELN QEL NQE
        |   0 113 114 128 129 227 242 242 257 355 356 370 371
    ----|----------------------------------------------------
    0   |
    113 | 113
    114 | 114   1
    128 | 128  15  14
    129 | 129  16  15   1
    227 | 227 114 113  99  98
    242 | 242 129 128 114 113  15
    242 | 242 129 128 114 113  15
    257 | 257 144 143 129 128  30  15  15
    355 | 355 242 241 227 226 128 113 113  98
    356 | 356 243 242 228 227 129 114 114  99  1
    370 | 370 257 256 242 241 143 128 128 113  15  14
    371 | 371 258 257 243 242 144 129 129 114  16  15   1
    484 | 484 371 370 356 355 257 242 242 227 129 128 114 113

        |  "" fal L   N   Q   LN  QE  fal LNQ ELN QEL NQE    fal: false
        |   0  99 113 114 128 227 257 299 355 356 370 371
    ----|------------------------------------------------
    0   |
    99  |  99
    113 | 113  14
    114 | 114  15   1
    128 | 128  29  15  14
    227 | 227 128 114 113  99
    257 | 257 158 144 143 129  30
    299 | 299 200 186 185 171  72  42
    355 | 355 256 242 241 227 128  98  56
    356 | 356 257 243 242 228 129  99  57   1
    370 | 370 271 257 256 242 143 113  71  15  14
    371 | 371 272 258 257 243 144 114  72  16  15   1
    484 | 484 385 371 370 356 257 227 185 129 128 114 113

    (Top) Spectral convolution for the theoretical spectrum of NQEL.
    The most frequent elements in the convolution between 57 and 200 are (multiplicities in parentheses):
    113 (8), 114 (8), 128 (8), 129 (8).
    (Bottom) Spectral convolution for the simulated spectrum of NQEL.
    The most frequent elements in the convolution between 57 and 200 are (multiplicities in parentheses):
    113 (4), 114 (4), 128 (4), 99 (3), 129 (3).

As predicted, some of the values in Figure 1 appear more frequently than others.
For example, 113 (the mass of "L") appears eight times in the convolution of the theoretical spectrum of "NQEL";
we say that 113 has multiplicity 8.
Six of the eight occurrences of 113 correspond to subpeptide pairs differing
in an "L": "L" and ""; "LN" and "N"; "EL" and "E"; "LNQ" and "NQ"; "QEL" and "QE"; "NQEL" and "NQE".

Spectral Convolution Problem
Compute the convolution of a spectrum.

Given:
A collection of integers Spectrum.

Return:
The list of elements in the convolution of Spectrum in decreasing order of their multiplicities.
If an element has multiplicity k, it should appear exactly k times.

Sample Dataset
0 137 186 323

Sample Output
137 137 186 186 323 49

═════════════════════════════════════════════════

Info.
  * sample dataset
        |   0 137 186 323
    ----|---------------
    0   |
    137 | 137
    186 | 186  49
    323 | 323 186 137      137(2), 186(2), 323(1), 49(1)

    137 137 186 186  323  49 --- a.  TODO: Notice that there's no zero!
    186 186 137 137  49  323 --- b.  TODO: Which is right?
    ---2--- ---2---  -1- -1-

  * top N elements with tie?
    [5,4,4,3,3,2,2,1,1]
    top 1 element with tie: 5
    top 2 element with tie: 5,4,4
    top 3 element with tie: 5,4,4,3,3

═════════════════════════════════════════════════

References:
*/
// #!/usr/bin/env javascript

// BA04H: Generate the Convolution of a Spectrum
function convolution(spectrum) {
    /**
     * [int] -> [int]
     * returns a convolution of spectrum in decreasing order
     */
    // sort unordered spectrum by ascending order
    spectrum.sort((a, b) => a - b);
    // get dictionary
    let freq_obj = {};
    for (let i = 1; i < spectrum.length; i++) {
        for (let j = 0; j < i; j++) {
            const diff = spectrum[i] - spectrum[j];
            if (!(diff in freq_obj)) freq_obj[diff] = 1;
            else freq_obj[diff]++;
        }
    }
    // convert dictionary to list for sorting
    let freq_list = [];
    for (const key in freq_obj) {
        freq_list.push([parseInt(key), freq_obj[key]]);
    }
    freq_list.sort((a, b) => b[1] - a[1]);
    // formatting
    let answer = []
    for (const tup of freq_list) {
        if (tup[0] === 0) continue;  // <-- Be careful! No zeros!
        answer.push(...Array(tup[1]).fill(tup[0]));
    }
    return answer;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04h.txt').toString().split("\n");
    const spectrum = lines[0].split(' ').map(Number);

    const startTime = performance.now();
    const conv = convolution(spectrum);
    console.log(conv.join(' '));
    console.log(`${performance.now() - startTime} milliseconds`);


    var stream = fs.createWriteStream("/home/wsl/rosalind/data/ba04h_output.txt");
    stream.once('open', function (fd) {
        stream.write(conv.join(' '));
        stream.end();
    });
}

// execute main function
main()