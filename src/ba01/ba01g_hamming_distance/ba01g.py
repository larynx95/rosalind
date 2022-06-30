"""
Rosalind: BA1G
Compute the Hamming Distance Between Two Strings

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
For example, CGAAT and CGGAC have two mismatches. The number of mismatches
between strings p and q is called the Hamming distance between these strings and
is denoted HammingDistance(p, q).

Hamming Distance Problem
Compute the Hamming distance between two DNA strings.

Given: Two DNA strings.

Return: An integer value representing the Hamming distance.

Sample Dataset
GGGCCGTTGGT
GGACCGTTGAC

Sample Output
3

═════════════════════════════════════════════════

    [ Where am I? ]

    * Where in the genome does DNA replication begin? (Chapter 01)
      -> Locating replication origin (OriC reion - DnaA box)
      -> minimum skewness (BA1F)
      -> Finding hidden message (DnaA box) in OriC
      ↓
    * Most frequent words problem
      -> count words (BA1A)
      -> find frequent words in a string (BA1B)
      -> find all occurrence of a pattern in a string (BA1D)
      -> PREV: Clump finidng problem (BA1E)
      ↓
    * Most frequent words with mismatch problem
      -> HERE: hamming distance (BA1G)
      -> NEXT: Find All Approximate Occurrences of a Pattern in a String (BA1H)

"""
# BA01G: hamming distance
def hdist(str1, str2):
    """
    (str,str) -> int
    """
    count = 0
    for i in range(len(str1.strip())):  # remove the last '\n'
        if str1[i] != str2[i]:
            count += 1
    return count


# BA01G: hamming distance
def hdist2(str1, str2):
    """
    (str,str) -> int
    """
    return sum([str1[i] != str2[i] for i in range(len(str1.strip()))])


# main function
def main():
    f = open('/home/wsl/rosalind/data/ba01g.txt', 'r')
    lines = f.readlines()
    str1 = lines[0]
    str2 = lines[1]
    f.close()

    print(hdist(str1, str2))
    print(hdist2(str1, str2))


# execute
if __name__ == "__main__":
    main()