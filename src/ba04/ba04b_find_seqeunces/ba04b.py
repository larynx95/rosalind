"""
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
- DNA is double helical structure.
  Are both of the two strands be transcribed as in this exercise?
  DNA Transcription (TODO: Is this exercise wrong? Yes it is.)
  https://www.nature.com/scitable/topicpage/dna-transcription-426/
  DNA is double-stranded, but only one strand serves as a template for transcription at any given time.
  This template strand is called the noncoding strand.
  The nontemplate strand is referred to as the coding strand
  because its sequence will be the same as that of the new RNA molecule.
- range, step
  https://www.programiz.com/python-programming/methods/built-in/range
- Remove sublist from list
  https://stackoverflow.com/questions/49737657/remove-sublist-from-list
- How do you split a list into evenly sized chunks?
  https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
- Print a list in reverse order with range()?
  https://stackoverflow.com/questions/7286365/print-a-list-in-reverse-order-with-range
- Coding strand
  https://en.wikipedia.org/wiki/Coding_strand
  coding strand (DNA):
    the DNA strand whose base sequence is identical to the base sequence of the RNA transcript produced
    (although with thymine replaced by uracil).
"""

from audioop import reverse
import time

ACIDs = {'F':['UUU','UUC'],'L':['UUA','UUG','CUU','CUC','CUA','CUG'],'S':['UCU','UCC','UCA','UCG','AGU','AGC'],'Y':['UAU','UAC'],'*':['UAA','UAG','UGA'],'C':['UGU','UGC'],'W':['UGG'],'P':['CCU','CCC','CCA','CCG'],'H':['CAU','CAC'],'Q':['CAA','CAG'],'R':['CGU','CGC','CGA','CGG','AGA','AGG'],'V':['GUU','GUC','GUA','GUG'],'A':['GCU','GCC','GCA','GCG'],'D':['GAU','GAC'],'E':['GAA','GAG'],'G':['GGU','GGC','GGA','GGG'],'I':['AUU','AUC','AUA'],'M':['AUG'],'T':['ACU','ACC','ACA','ACG'],'N':['AAU','AAC'],'K':['AAA','AAG']}
CODONs = {'UUU':'F','CUU':'L','AUU':'I','GUU':'V','UUC':'F','CUC':'L','AUC':'I','GUC':'V','UUA':'L','CUA':'L','AUA':'I','GUA':'V','UUG':'L','CUG':'L','AUG':'M','GUG':'V','UCU':'S','CCU':'P','ACU':'T','GCU':'A','UCC':'S','CCC':'P','ACC':'T','GCC':'A','UCA':'S','CCA':'P','ACA':'T','GCA':'A','UCG':'S','CCG':'P','ACG':'T','GCG':'A','UAU':'Y','CAU':'H','AAU':'N','GAU':'D','UAC':'Y','CAC':'H','UAG':'stop','GAC':'D','CAA':'Q','AAA':'K','GAA':'E','AAC':'N','CAG':'Q','AAG':'K','GAG':'E','UGU':'C','UAA':'stop','CGU':'R','AGU':'S','GGU':'G','UGC':'C','CGC':'R','AGC':'S','GGC':'G','CGA':'R','AGA':'R','UGA':'stop','GGA':'G','UGG':'W','CGG':'R','AGG':'R','GGG':'G'}


def reverse_complement(dna):
    """
    reverse complement
    >>> reverse_complement('ACCGT')
        ACGGT
    """
    comp = ''
    for nuc in dna[::-1]:
        if nuc == 'A':
            comp += 'T'
        if nuc == 'C':
            comp += 'G'
        if nuc == 'G':
            comp += 'C'
        if nuc == 'T':
            comp += 'A'
    return comp


def transcribe(DNA):
    """
    DNA transcription
    >>> transcribe('ACGTTC')
        ACTUUC
    """
    mRNA = ''
    for nuc in DNA:
        if nuc == 'T':
            mRNA += 'U'
        else:
            mRNA += nuc
    return mRNA


def rev_transcribe(mRNA):
    DNA = ''
    for rn in mRNA:
        if rn == 'U':
            DNA += 'T'
        else:
            DNA += rn
    return DNA


def translate(mRNA):
    """
    translate mRNA to protein (amino acids)
    Suppose that the direction of mRNA is 5 to 3.
    >>> translate('ACGUAC)
        TY
    """
    protein = ''
    for i in range(0, len(mRNA), 3):
        if len(mRNA) < i+3:
            return protein
        codon = mRNA[i:i+3]
        aa = CODONs[codon]
        protein += aa
    return protein

def find_subs_WRONG(DNA, AA):
    """
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    aa = 'MA'
    >>> find_subs_WRONG(dna, aa)
        ['ATGGCC', 'ATGGCC']
    """
    result = []
    mRNAs = [transcribe(DNA), transcribe(reverse_complement(DNA))]
    for mRNA in mRNAs:
        idx_last = ((len(mRNA) - 3 * len(AA)) // 3) * 3
        for i in range(0, idx_last + 1, 3 * len(AA)):  # <--- this is wrong
            codons = mRNA[i:i+3*len(AA)]
            aacids = translate(codons)
            if aacids == AA:
                result.append(rev_transcribe(codons))
    return result

def find_subs_WRONG2(DNA, AA):
    """
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    aa = 'MA'
    >>> find_subs_WRONG2(dna, aa)
        ['ATGGCC', 'ATGGCC', 'ATGGCC']
    """
    result = []
    mRNAs = [transcribe(DNA), transcribe(reverse_complement(DNA))]
    for mRNA in mRNAs:  # <-- this is wrong
        idx_last = ((len(mRNA) - 3 * len(AA)) // 3) * 3
        for i in range(0, idx_last + 1):
            codons = mRNA[i:i+3*len(AA)]
            aacids = translate(codons)
            if aacids == AA:
                result.append(rev_transcribe(codons))  # <--- this is wrong
    return result

def find_subs(DNA, AA):
    """
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    aa = 'MA'
    >>> find_subs_WRONG(dna, aa)
        ['ATGGCC', 'ATGGCC', 'GGCCAT']
    """
    result = []
    rev = False
    mRNAs = [transcribe(DNA), transcribe(reverse_complement(DNA))]
    for j in range(len(mRNAs)):
        if j == 1:
            rev = True
        idx_last = ((len(mRNAs[j]) - 3 * len(AA)) // 3) * 3
        for i in range(0, idx_last + 1):
            codons = mRNAs[j][i:i+3*len(AA)]
            aacids = translate(codons)
            if aacids == AA:
                seq_DNA = rev_transcribe(codons)
                if rev:
                    seq_DNA = reverse_complement(seq_DNA)
                result.append(seq_DNA)
    return result

def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.
    https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main():
    f = open('/home/wsl/rosalind/data/ba04b.txt', 'r')
    DNA, AA = [line.strip() for line in f.readlines()]
    f.close()

    print(reverse_complement('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'))
    start_time = time.time()
    """
    subs = find_subs(DNA, AA)
    for sub in subs:
        print(sub)
    """
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()
