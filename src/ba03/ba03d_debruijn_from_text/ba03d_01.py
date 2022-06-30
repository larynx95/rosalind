'''
Rosalind: BA3D
Construct the "De Bruijn[broin]" Graph of a String

Given a genome Text, PathGraph_k(Text) is the path consisting of |Text| - k + 1 edges,
where the i-th edge of this path is labeled by the i-th k-mer in Text
and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text.
The de Bruijn graph DeBruijn_k(Text) is formed
by gluing identically labeled nodes in PathGraph_k(Text).

De Bruijn Graph from a String Problem
Construct the de Bruijn graph of a string.

Given: An integer k and a string Text.

Return:DeBruijnk(Text), in the form of an adjacency list.

Sample Dataset
4
AAGATTCTCTAC

Sample Output
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC
TTC -> TCT

-------------------------------------------------

solution by apense
https://rosalind.info/users/apense/
'''

import sys
import time

def kmer_generator(dna, kmer_size):
    position = 0
    while(position+kmer_size <= len(dna)):
        yield dna[position:position+kmer_size]
        position += 1

def de_bruijn(dna, size):
    result = {}
    kmers = [m for m in kmer_generator(dna, size-1)]
    for i in range(len(kmers)-1):
        try:
            result[kmers[i]].append(kmers[i+1])
        except KeyError:
            result[kmers[i]] = [kmers[i+1]]
    return result

def main():
    f = open('/home/wsl/rosalind/data/ba03d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    text = lines[1]
    f.close()

    start_time = time.time()
    m = de_bruijn(text, k)
    original_stdout = sys.stdout
    with open('/home/wsl/rosalind/data/ba03d_output.txt', 'w') as f:
        sys.stdout = f
        for key in m:
            f.write(key + " -> " + ','.join([v for v in m[key]]) + "\n")
        sys.stdout = original_stdout
    print("--- %s seconds ---" % (time.time() - start_time))

# main function
if __name__ == "__main__":
    main()
