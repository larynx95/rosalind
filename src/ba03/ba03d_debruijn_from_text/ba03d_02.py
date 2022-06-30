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

solution by Md. Ahsan Kabir
https://rosalind.info/users/shohagh_sust/
'''

import sys
import time

def Kmer(k,Text):
     s = set()
     for i in range(0,len(Text)-k+1):
         s.add(Text[i:i+k])
     return s

def posFind(s,Text):
    pos = list()
    for i in range(0,len(Text)-len(s)):
        if s==Text[i:i+len(s)]:
            pos.append(i)
    return pos

def DeBruijnk(k,Text):
    kmer = Kmer(k-1,Text)
    for i in kmer:
        cur = i
        pos = posFind(cur,Text)
        s = ""
        for j in range(0,min(1,len(pos))):
            s+=" "+Text[pos[j]+1:pos[j]+1+k-1]
        for j in range(1,len(pos)):
            s += ","+Text[pos[j]+1:pos[j]+1+k-1];
        if s!="":
            s = cur+" ->"+s
            print(s)
    return

def main():
    f = open('/home/wsl/rosalind/data/ba03d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    text = lines[1]
    f.close()

    start_time = time.time()
    original_stdout = sys.stdout
    with open('/home/wsl/rosalind/data/ba03d_output.txt', 'w') as f:
        sys.stdout = f
        DeBruijnk(k,text)
        sys.stdout = original_stdout
    print("--- %s seconds ---" % (time.time() - start_time))


# main function
if __name__ == "__main__":
    main()