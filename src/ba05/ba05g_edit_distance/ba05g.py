"""
Rosalind: BA5G
Compute the Edit Distance Between Two Strings

In 1966, Vladimir Levenshtein introduced the notion of the edit distance between two strings
as the minimum number of edit operations needed to transform one string into another.
Here, an edit operation is the insertion, deletion, or substitution of a single symbol.
For example, TGCATAT can be transformed into ATCCGAT with five edit operations,
implying that the edit distance between these strings is at most 5.

Edit Distance Problem
Find the edit distance between two strings.

Given: Two amino acid strings.

Return: The edit distance between these strings.

Sample Dataset
PLEASANTLY
MEANLY

Sample Output
5

═════════════════════════════════════════════════

Plan 1.
  * sample dataset
    PLEASANTLY
    -MEA--N-LY
    12  34 5
  * This exercise is similar to the longest common sequence (lcs, BA5C->BA5E) problem

          M   E   A   N   L   Y   ↘: match or mismatch
        ┌───┬───┬───┬───┬───┬───┐ → : insertion
      P ↓   │   │   │   │   │   │ ↓ : deletion
        ├───┼───┼───┼───┼───┼───┤
      L │↘ │   │   │   │   │   │ modification:
        ├───┼───┼───┼───┼───┼───┤     mismatch : +1 ↘
      E │   │↘ │   │   │   │   │     deletion : +1 ↓
        ├───┼───┼───┼───┼───┼───┤     insertion: +1 →
      A │   │   │↘ │   │   │   │
        ├───┼───┼───┼───┼───┼───┤
      S │   │   │   ↓   │   │   │
        ├───┼───┼───┼───┼───┼───┤
      A │   │   │   ↓   │   │   │
        ├───┼───┼───┼───┼───┼───┤
      N │   │   │   │↘ │   │   │
        ├───┼───┼───┼───┼───┼───┤
      T │   │   │   │   ↓   │   │
        ├───┼───┼───┼───┼───┼───┤
      L │   │   │   │   │↘ │   │
        ├───┼───┼───┼───┼───┼───┤
      Y │   │   │   │   │   │↘ │
        └───┴───┴───┴───┴───┴───┘

═════════════════════════════════════════════════

References:
"""
#!/usr/bin/env python
import time


def edit_distance(v, w):
    """
    (a,a) -> int
    returns a set of backtrack
    >>> edit_distance(edges_blosum62, 'PLEASANTLY','MEANLY')
        8
    """
    s = {}
    for i in range(len(v)+1):
        s[(i,0)] = i
    for j in range(len(w)+1):
        s[(0,j)] = j
    for i in range(len(v)):
        for j in range(len(w)):
            insertion = s[(i,j+1)] + 1                         # insertion point
            deletion  = s[(i+1,j)] + 1                         # deletion point
            mismatch  = s[(i,j)] + (1 if v[i] != w[j] else 0)  # mismatch point
            s[(i+1,j+1)] = min(insertion, deletion, mismatch)
    return s[(len(v),len(w))]


def main():
    f = open('/home/wsl/rosalind/data/ba05g.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    astr = lines[0]
    bstr = lines[1]
    f.close()

    start_time = time.time()
    print(edit_distance(astr, bstr))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
