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

solution by Matthew Burch
https://rosalind.info/users/cryptc/
'''

import sys
import time

def DeBruijn_k(k, dna):
    kmerLength = k-1
    graph = {}
    for i in range(1,len(dna)-k+1+1):
        prevKmer = dna[i-1:i-1+kmerLength]
        kmer = dna[i:i+kmerLength]
        graph[prevKmer] = graph.get(prevKmer, [])
        if prevKmer[1:] == kmer[:-1]:
            graph[prevKmer].append(kmer)
    return graph

def main():
    f = open('/home/wsl/rosalind/data/ba03d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    text = lines[1]
    f.close()

    graph = DeBruijn_k(k,text)
    print(graph)  # dictionary

    '''
    start_time = time.time()
    original_stdout = sys.stdout
    with open('../data/ba3d_output.txt', 'w') as f:
        sys.stdout = f
        pass
        sys.stdout = original_stdout
    print("--- %s seconds ---" % (time.time() - start_time))
    '''

# main function
if __name__ == "__main__":
    main()