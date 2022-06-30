/*
Rosalind: BA4C
Generate the Theoretical Spectrum of a Cyclic Peptide

The workhorse of peptide sequencing is the mass spectrometer,
an expensive molecular scale that shatters molecules into pieces and then weighs the resulting fragments.
The mass spectrometer measures the mass of a molecule in daltons (Da);
1 Da is approximately equal to the mass of a single nuclear particle (i.e., a proton or neutron).

We will approximate the mass of a molecule by simply adding the number of protons
and neutrons found in the molecule’s constituent atoms,
which yields the molecule’s integer mass.
For example, the amino acid "Gly", which has chemical formula C2H3ON,
has an integer mass of 57, since 2·12 + 3·1 + 1·16 + 1·14 = 57.
Yet 1 Da is not exactly equal to the mass of a proton/neutron,
and we may need to account for different naturally occurring isotopes of each atom when weighing a molecule.
 As a result, amino acids typically have non-integer masses
 (e.g., "Gly" has total mass equal to approximately 57.02 Da);
 for simplicity, however, we will work with the integer mass table given in Figure 1.

  G  A  S  P  V  T   C   I   L   N   D   K   Q   E   M   H   F   R   Y   W
  57 71 87 97 99 101 103 113 113 114 115 128 128 129 131 137 147 156 163 186

The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide),
is the collection of all of the masses of its subpeptides,
in addition to the mass 0 and the mass of the entire peptide.
We will assume that the theoretical spectrum can contain duplicate elements,
as is the case for "NQEL" (shown in Figure 2), where "NQ" and "EL" have the same mass.

    L   N   Q   E   LN  NQ  EL  QE  LNQ ELN QEL NQE NQEL  <-- cyclic
  0 113 114 128 129 227 242 242 257 355 356 370 371 484

Generating Theoretical Spectrum Problem
Generate the theoretical spectrum of a cyclic peptide.

Given: An amino acid string Peptide.

Return: Cyclospectrum(Peptide).

Sample Dataset
LEQN

Sample Output
0 113 114 128 129 227 242 242 257 355 356 370 371 484

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
      -> LinearSpectrum (BA4J)
      -> LinearScore (BA4K)
      -> PREV: the number of Linear peptides (BA4D)
      ↓
    * Cyclo
      -> HERE: Cyclospectrum (BA4C)
      -> NEXT: CycloScore (BA4F)

    [ summary ]

    # Linear vs. Cyclo

    Linear-peptides, spectrum ┌ Peptides -> Spectrum
                              │   ├ Accumulating
                              │   └ Reusing (or Referencing)
                              └ Masses -> Spectrum
    Cyclo-peptides, spectrum  ┌ Peptides -> Spectrum
                              │   ├ Accumulating
                              │   └ Reusing (or Referencing)
                              └ Masses -> Spectrum

    # Chain of fuctions

    linear_subpeptides -> spectrum_linear_from_subpeptides
    cyclo-peptides -> spectrum_cyclo_from_subpeptides


Plan 1.
  * peptides, then spectrum

    1. Linear-subpeptides, spectrum
      1) peptides, then spectrum
      - constructing list of subpeptides

        given peptide: NQEL
        i  j  append
        ----------------------- --> generalized rule
        0  0  peptide[0:1] N        outer loop: i in range(len(peptide))
           1  peptide[1:2] Q        inner loop: j in range(len(peptide)-i)
           2  peptide[2:3] E
           3  peptide[3:4] L        subpeptide = peptide[j:i+j+1]
        1  0  peptide[0:2] NQ
           1  peptide[1:3] QE
           2  peptide[2:4] EL
        2  0  peptide[0:3] NQE
           1  peptide[1:4] QEL
        3  0  peptide[0:4] NQEL

    - constructing list of subpeptides, masses

        N   Q   E   L   NQ  QE  EL  NQE QEL NQEL
      0 114 128 129 113 242 257 242 371 370 484
      after sorting
      0 113 114 128 129 242 242 257 370 371 484

      Mass(QE) = 128 + 129 = 257
      PrefixMass = 0,114,242,371,484
                   values in above list were calculated as belows...
                   0,N,  NQ, NQE,NQEL
                   0, N        , NQ             , NQE                    , NQEL
                   0, 0+Mass(N), Mass(N)+Mass(Q), Mass(N)+Mass(Q)+Mass(E), Mass(N)+Mass(Q)+Mass(E)+Mass(L)
                   0, 0+114    , 114+128=242    , 242+129=371            , 371+113=484
                   0, 114      , 242            , 371                    , 484

      i  peptides  peptide[i:j]             Mass(peptide)
        --------------------------------------------------------------------------
      0  N Q E L   [0:1] [1:2] [2:3] [3:4]  p[1]-p[0] p[2]-p[1] p[3]-p[2] p[4]-p[3]
      1  NQ QE EL  [0:2] [1:3] [2:4]                  p[2]-p[0] p[3]-p[1] p[4]-p[2]
      2  NQE QEL   [0:3] [1:4]                                  p[3]-p[0] p[4]-p[1]
      3  NQEL      [0:4]                                                  p[4]-p[0]

            subpeptides = []                     spectrum = [0]
            for i=0 to len(p):                   for i=0 for len(p):
                for j=i+1 to len(p)+1:               for j=0 for len(p)-i:
                    subpeptides.append(p[0:j])           spectrum.append(Pref[j+i+1]-Pref[j])

      2) masses, then spectrum

    2. cyclo-subpeptides, spectrum
      1) peptides, then spectrum
      (1) constructing cyclo-spectrum list of subpeptides - "accumutating" a result list

      focusing on the letter of amino acid
        NQEL  ''                                   common operation (updating result list)
              N Q E L          '' + N,Q,E,L        add the next aacid to each element
              NQ QE EL LN      N+Q Q+E E+L L+N     add the next aacid to each element
              NQE QEL ELN LNQ  NQ+E QE+L EL+N LN+Q add the next aacid to each element
              NQEL             NQE+L               add the next aacid to each element

      pt='NQEL...', len(pt) == n, ls=['',N,Q,E,L,...]  <-- list of peptides with just one amino acid
      i=0 ['',N,Q,E,L]
          ls[-4]+pt[0], ls[-3]+pt[1], ls[-2]+pt[2], ls[-1]+pt[3]
      i=1 ['',N,Q,E,L,NQ,QE,EL,LN]
          ls[-4]+pt[1], ls[-3]+pt[2], ls[-2]+pt[3], ls[-1]+pt[0]
      i=2 ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ]
          ls[-4]+pt[2], ls[-3]+pt[3], ls[-2]+pt[0], ls[-1]+pt[1]
      i=k ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ, ...]
          ls[-n]+pt[(k+0)%n], ls[-n+1]+pt[(k+1)%n], ... , ls[-n+(n-2)]+pt[(k+(n-2))%n] , ls[-n+(n-1)]+pt[(k+(n-1))%n]
      i=n ['',N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ, ...]
          ls+pt

      (2) constructing cyclo-spectrum list of subpeptides - "reusing" a string repeatedly ("NQEL")

        L   N   Q   E   LN  NQ  EL  QE  LNQ ELN QEL NQE NQEL  <-- cyclic
      0 113 114 128 129 227 242 242 257 355 356 370 371 484

      i  peptides         peptide[i:j]       generalized rule
      --------------------------------------------------------------------------
      0  N Q E L          [0]                if i+j+1 <= 4, [j:i+j+1]
                          [1]                if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                          [2]
                          [3]
      1  NQ QE EL LN      [0:2]              if i+j+1 <= 4, [j:i+j+1]
                          [1:3]              if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                          [2:4]
                          [3:]    +[:1]
      2  NQE QEL ELN LNQ  [0:3]              if i+j+1 <= 4, [j:i+j+1]
                          [1:4]              if i+j+1 > 4,  [j:] + [(i+j+1)%4]
                          [2:]    +[:1]
                          [3:]    +[:2]
      3  NQEL             [0:4]

      2) masses, then spectrum

Plan 2.
  * textbook alrogorithms - charging station (textbook pp.211~212)
  * creating a "PrefixMass" list, then get spectrum
    PrefixMass('NQEL') = [0,114,242,371,484]
                            N   NQ  NQE NQEL
  1. linear
  - Pseudocode:
    ╔═══════════════════════════════════════════════════════════════════════╗
    ║ LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)                     ║
    ║     # A. section for constructing PrefixMass list                     ║
    ║     PrefixMass(0) <- 0                                                ║
    ║     for i <- 1 to |Peptide|                                           ║
    ║         for j <- 1 to 20                                              ║
    ║             if AminoAcid(j) =  i-th amino acid in Peptide             ║
    ║                 PrefixMass(i) <- PrefixMass(i - 1) + AminoAcidMass(j) ║
    ║     # B. section for constructing spectrum                            ║
    ║     LinearSpectrum <- a list consisting of the single integer 0       ║
    ║     for i <- 0 to |Peptide| - 1                                       ║
    ║         for j <- i + 1 to |Peptide|                                   ║
    ║             add PrefixMass(j) - PrefixMass(i) to LinearSpectrum       ║
    ║     return sorted list LinearSpectrum                                 ║
    ╚═══════════════════════════════════════════════════════════════════════╝

    i j   [0,113,242,370,484]
    -------------------------
    0 1-4  0,113,242,370,484
    1 2-4        128,257,371
    2 3-4            128,242
    3 4                  114

  2. cyclo
  - pseudocode (textbook):
    ╔═════════════════════════════════════════════════════════════════════════════════════╗
    ║ CYCLICSPECTRUM(Peptide, AMINOACID, AMINOACIDMASS)                                   ║
    ║     # A. section for constructing PrefixMass list                                   ║
    ║     PREFIXMASS(0) 0                                                                 ║
    ║     for i <- 1 to Peptide                                                           ║
    ║         for j <- 1 to 20                                                            ║
    ║             if AMINOACID(j) = i-th amino acid in Peptide                            ║
    ║                 PREFIXMASS(i) PREFIXMASS(i - 1) + AMINOACIDMASS(j)                  ║
    ║     peptideMass <- PREFIXMASS(|Peptide|)                                            ║
    ║     # B. section for constructing spectrum list                                     ║
    ║     CyclicSpectrum <- a list consisting of the single integer 0                     ║
    ║     for i <- 0 to |Peptide| - 1                                                     ║
    ║         for j <- i + 1 to |Peptide|                                                 ║
    ║             add PREFIXMASS(j) - PREFIXMASS(i) to CyclicSpectrum                     ║
    ║             if i > 0 and j < |Peptide|                                              ║
    ║                 add peptideMass - (PREFIXMASS(j) - PREFIXMASS(i)) to CyclicSpectrum ║
    ║     return sorted list CyclicSpectrum                                               ║
    ╚═════════════════════════════════════════════════════════════════════════════════════╝

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
    return lspec.sort((a, b) => a - b);  // <-- 'sort()' sorts only string. Be careful.
}

// BA04C: cyclic (cyclo) spectrum
function spectrum_cyclo(peptide) {
    /**
     * str -> [int]
     * returns a cyclospectrum
     */
    let prefixmass = [0];
    for (const acid of peptide) {
        prefixmass.push(prefixmass[prefixmass.length - 1] + weights[acid]);
    }
    const mass = prefixmass[prefixmass.length - 1];
    let cspec = [0];
    for (let i = 0; i < peptide.length; i++) {
        for (let j = i + 1; j < peptide.length + 1; j++) {
            cspec.push(prefixmass[j] - prefixmass[i]);
            if (i > 0 && j < peptide.length) {
                cspec.push(mass - (prefixmass[j] - prefixmass[i]));
            }
        }
    }
    return cspec.sort((a, b) => a - b);  // <-- 'sort()' sorts only string. Be careful.
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04c.txt').toString().split("\n");
    const peptide = lines[0];

    const startTime = performance.now();
    const cspec = spectrum_cyclo(peptide);
    console.log(cspec.join(' '));
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()