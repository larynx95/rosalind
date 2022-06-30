'''
Rosalind: BA4E
Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

In "Compute the Number of Peptides of Given Total Mass",
we first encountered the problem of reconstructing a cyclic peptide from its theoretical spectrum;
this problem is called the Cyclopeptide Sequencing Problem and is given below.
It is solved by the following algorithm.

    CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides <- a set containing only the empty peptide
        while Peptides is nonempty
            Peptides <- Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides

Cyclopeptide Sequencing Problem
Given an ideal experimental spectrum,
find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.

Return: Every amino acid string Peptide such that
        Cyclospectrum(Peptide) = Spectrum (if such a string exists).

Sample Dataset
0 113 128 186 241 299 314 427

Sample Output
186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

-------------------------------------------------

solution by NathanielLovin
https://github.com/NathanielLovin/Rosalind/blob/master/BA4E.py
'''

def cyclospectrum(peptide):
	spec = [0]
	for x in range(1,len(peptide)):
		for i in range(len(peptide)):
			if i+x >= len(peptide):
				y = i+x-len(peptide)
				spec.append(sum(peptide[i:])+sum(peptide[:y]))
			else:
				spec.append(sum(peptide[i:i+x]))
	spec.append(sum(peptide))
	spec.sort()
	return spec


def linospectrum(peptide):
	spec = [0]
	for x in range(1,len(peptide)):
		for i in range(len(peptide)):
			if i+x <= len(peptide):
				spec.append(sum(peptide[i:i+x]))
	spec.append(sum(peptide))
	spec.sort()
	return spec


def expand(ls):
	peps = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
	newlist = []
	for x in ls:
		for z in peps:
			y = x.copy()
			y.append(z)
			newlist.append(y)
	return newlist


def output(pep):
	p = ""
	for x in pep:
		p += str(x) + '-'
	print(p[:-1])


def cyclopeptideSequencing(spectrum):
	peptides = [[]]
	while peptides != []:
		peptides = expand(peptides)
		newpeptides = []
		for pep in peptides:
			if sum(pep) == spectrum[-1]:
				spec = cyclospectrum(pep)
				if len(spec) == len(spectrum):
					match = True
					for i in range(len(spec)):
						if spec[i] != spectrum[i]:
							match = False
							break;
					if match:
						output(pep)
			else:
				spec = linospectrum(pep)
				match = True
				i = -1
				for x in spec:
					xmatch = False
					for y in range(i+1,len(spectrum)):
						if x==spectrum[y]:
							xmatch = True
							i = y
							break;
					if xmatch == False:
						match = False
						break;
				if match:
					newpeptides.append(pep)
		peptides = newpeptides

spectrum = [0, 113, 115, 115, 115, 129, 131, 156, 186, 186, 230, 242, 244, 246, 271, 299, 301, 317, 342, 357, 359, 361, 414, 428, 432, 457, 457, 472, 473, 490, 543, 543, 547, 570, 588, 588, 603, 643, 658, 658, 676, 699, 703, 703, 756, 773, 774, 789, 789, 814, 818, 832, 885, 887, 889, 904, 929, 945, 947, 975, 1000, 1002, 1004, 1016, 1060, 1060, 1090, 1115, 1117, 1131, 1131, 1131, 1133, 1246]
cyclopeptideSequencing(spectrum)