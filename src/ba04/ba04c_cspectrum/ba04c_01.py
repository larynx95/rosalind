'''
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

-------------------------------------------------

by NathanielLovin
https://github.com/NathanielLovin/Rosalind/blob/master/BA4C.py
'''

mass = {}
mass['G'] = 57
mass['A'] = 71
mass['S'] = 87
mass['P'] = 97
mass['V'] = 99
mass['T'] = 101
mass['C'] = 103
mass['I'] = 113
mass['L'] = 113
mass['N'] = 114
mass['D'] = 115
mass['K'] = 128
mass['Q'] = 128
mass['E'] = 129
mass['M'] = 131
mass['H'] = 137
mass['F'] = 147
mass['R'] = 156
mass['Y'] = 163
mass['W'] = 186

def proteinMass(pep):
	weight = 0
	for x in pep:
		weight += mass[x]
	return weight

peptide = "SVYQERFHWNIHSH"
spec = [0]
for x in range(1,len(peptide)):
	for i in range(len(peptide)):
		if i+x > len(peptide):
			y = i+x-len(peptide)
			spec.append(proteinMass(peptide[i:]+peptide[:y]))
		else:
			spec.append(proteinMass(peptide[i:i+x]))
spec.append(proteinMass(peptide))
spec.sort()
p = ""
for x in spec:
	p +=  str(x) + " "
print(p)