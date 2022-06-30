"""
Rosalind: BA5C
Find a Longest Common Subsequence of Two Strings

Longest Common Subsequence Problem
Given: Two strings.

Return: A longest common subsequence of these strings.

Sample Dataset
AACCTTGG
ACACTGTGA

Sample Output
AACTGG

═════════════════════════════════════════════════

solution by Carlos Garcia
"""

f = open('/home/wsl/rosalind/data/ba05c.txt', 'r')
a = f.read()
a = a.split('\n')
V = a[0]
W = a[1]


def LCSBacktrack(v,w):
    s = {}
    Backtrack = {}
    for i in range(len(v)+1):
        s[(i,0)]= 0
    for j in range(len(w)+1):
        s[(0,j)] = 0
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                s[(i+1,j+1)] = s[(i,j)] + 1
            else:
                candidates = [s[(i+1,j)],s[(i,j+1)]]
                s[(i+1,j+1)] = max(candidates)
    return s


BACKTRACK = LCSBacktrack(V,W)


def OUTPUTLcs(Backtrack,v,w):
    i = len(v)
    j = len(w)
    result = []
    while i*j !=0:
        if Backtrack[(i,j)] == Backtrack[(i-1,j)]:
            i -= 1
        elif Backtrack[(i,j)] == Backtrack[(i,j-1)]:
            j -= 1
        else:
            result.append(v[i-1])
            j -= 1
            i -= 1
    return ''.join(result[::-1])


print (OUTPUTLcs(BACKTRACK,V,W))
