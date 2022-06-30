/*
Rosalind: BA4B
Find Substrings of a Genome Encoding a Given Amino Acid String

There are three different ways to divide a DNA string into codons for translation,
one starting at each of the first three starting positions of the string.
These different ways of dividing a DNA string into codons are called reading frames.
Since DNA is double-stranded, a genome has six reading frames (three on each strand),
as shown in Figure 1.

We say that a DNA string Pattern encodes an amino acid string Peptide
if the RNA string transcribed from either Pattern
or its "reverse complement Pattern" translates into Peptide.

Peptide Encoding Problem
Find substrings of a genome encoding a given amino acid sequence.

Given: A DNA string Text and an amino acid string Peptide.

Return: All substrings of Text encoding Peptide (if any such substrings exist).

Sample Dataset
ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
MA

Sample Output
ATGGCC
GGCCAT
ATGGCC

═════════════════════════════════════════════════

    [ Where am I? ] How do we sequence antibiotics?

    * Central Dogma
      -> PREV: translation (BA4A)
      -> HERE: reverse translation (BA4B)
      ↓
    * Dodging the Central Dogma: Non-ribosomal peptides (NRPs)
      totally different from the Central Dogma
      ↓
    * Linear
      -> NEXT: LinearSpectrum (BA4J)

Plan 1.
- focusing on two DNA strings and their two mRNA strings
  DNA:  ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
  rDNA: TCACCCGTTAATACGGGTACTATTGATCTCAGTTCTGGGGGCCATGGCCAT
  TODO: Which DNA string is the "coding" strand?
  In this exercise, it is assumed that ...
    the given DNA string is a coding strand (5->3)
    and also assumed that the reverse complement DNA string is a coding strand. (???)
    (TODO: Is this assumption correct?)
- anyway...the next steps are...
  - get two DNA strings: DNA and its reverse complement string
  - get two mRNA strings, one from DNA and the other from its reverse complement
  - do transcription and translation and compare result with amino acids already given in this exercise

Plan 2.
- focusing only on DNA substrings
  amino acids: MA  --> length of codons: len(amino acids) * 3 = 6
  DNA: ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
       ------  here, create complement DNA, transcribe, translate, compare ...
        ------  repeat
         ------  repeat again ...
- TODO: Implement this plan.

═════════════════════════════════════════════════

References:
- central dogma
- DNA is double helical structure.
  Are both of the two strands be transcribed as in this exercise?
  DNA Transcription (TODO: Is this exercise wrong? Yes it is.)
  https://www.nature.com/scitable/topicpage/dna-transcription-426/
  DNA is double-stranded, but only one strand serves as a template for transcription at any given time.
  This template strand is called the noncoding strand.
  The nontemplate strand is referred to as the coding strand
  because its sequence will be the same as that of the new RNA molecule.
- Coding strand
  https://en.wikipedia.org/wiki/Coding_strand
  coding strand (DNA):
    the DNA strand whose base sequence is identical to the base sequence of the RNA transcript produced
    (although with thymine replaced by uracil).
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

// BA01C: reverse complement
function reverse_complement(dna) {
    /**
     * str -> str
     * returns a reverse complement (BA1C)
     * reverse_complement("ACGGTT") == AACCGT
     */
    dic = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' };
    dna = [...dna].reverse().join("");
    rev = "";
    for (const nuc of dna) {
        rev += dic[nuc];
    }
    return rev;
}

// BA04B: transcription
function transcribe(dna) {
    /**
     * str -> str
     * returns a mRNA from a DNA string
     * >>> transcribe("ACGT")
     *     ACGU
     */
    return Array.from(dna).map(nuc => nuc === 'T' ? 'U' : nuc).join('');
}

// BA04A: translation
function translate(mrna) {
    /**
     * str -> str
     * returns a amino acids from a mrna
     */
    let aa = "";
    for (let i = 0; i < mrna.length; i += 3) {
        const codon = mrna.slice(i, i + 3);
        if (['UAA', 'UAG', 'UGA'].includes(codon)) break;
        aa += codons[codon];
    }
    return aa;
}

// BA04B: find substrings
function find_substrs(dna, pep) {
    /**
     * (str,str) -> [str]
     * returns a list of substrings
     */
    let result = [];
    for (let i = 0; i < dna.length - 5; i++) {
        const subs = dna.slice(i, i + 6);
        const subs_rev = reverse_complement(subs);
        const mrna = transcribe(subs);
        const mrna_rev = transcribe(subs_rev);
        const aa = translate(mrna);
        const aa_rev = translate(mrna_rev);
        if (pep === aa || pep == aa_rev) result.push(subs);
    }
    return result;
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04b.txt').toString().split("\n");
    const dna = lines[0];
    const pep = lines[1];
    const startTime = performance.now();
    const substrs = find_substrs(dna, pep);
    for (const substr of substrs) {
        console.log(substr);
    }
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()