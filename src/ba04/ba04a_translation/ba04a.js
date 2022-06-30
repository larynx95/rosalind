/*
Rosalind: BA4A
Translate an RNA String into an Amino Acid String

Much like replication,
the chemical machinery underlying transcription and translation is fascinating,
but from a computational perspective, both processes are straightforward.
Transcription simply transforms a DNA string into an RNA string
by replacing all occurrences of "T" with "U".
The resulting strand of RNA is translated into an amino acid sequence via the genetic code;
this process converts each 3-mer of RNA, called a codon, into one of 20 amino acids.

As illustrated in Figure 1,
each of the 64 RNA codons encodes its own amino acid (some codons encode the same amino acid),
with the exception of three stop codons that do not translate into amino acids and serve to halt translation.
For example, the DNA string "TATACGAAA" transcribes into the RNA string "UAUACGAAA",
which in turn translates into the amino acid string "Tyr-Thr-Lys".

The following problem asks you to find the translation of an RNA string into an amino acid string.

Protein Translation Problem
Translate an RNA string into an amino acid string.

Given: An RNA string Pattern.

Return: The translation of Pattern into an amino acid string Peptide.

Sample Dataset
AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output
MAMAPRTEINSTRING

═════════════════════════════════════════════════

    [ Where am I? ]

    * Central Dogma
      -> HERE: translation (BA4A)
      -> NEXT: reverse translation (BA4B)

═════════════════════════════════════════════════

References:
-
*/
// #!/usr/bin/env javascript

const aacids = {
    'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Y': ['UAU', 'UAC'],
    '*': ['UAA', 'UAG', 'UGA'], 'C': ['UGU', 'UGC'], 'W': ['UGG'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'], 'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG']
}
const codons = {
    'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I',
    'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L',
    'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S',
    'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
    'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N',
    'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'UAG': 'stop', 'GAC': 'D', 'CAA': 'Q', 'AAA': 'K',
    'GAA': 'E', 'AAC': 'N', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'UGU': 'C', 'UAA': 'stop',
    'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
    'CGA': 'R', 'AGA': 'R', 'UGA': 'stop', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

// BA04A: translate
function translate(seq) {
    /**
     * str -> str
     * returns a peptide
     */
    let peptide = "";
    for (let i = 0; i < seq.length; i += 3) {
        const codon = seq.slice(i, i + 3);
        if (['UAA', 'UAG', 'UGA'].includes(codon)) break;
        peptide += codons[codon];
    }
    return peptide;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04a.txt').toString().split("\n");
    const seq = lines[0];

    const startTime = performance.now();
    console.log(translate(seq));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()