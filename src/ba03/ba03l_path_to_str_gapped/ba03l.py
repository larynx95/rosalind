"""
Rosalind: BA3L
Construct a String Spelled by a Gapped Genome Path

Gapped Genome Path String Problem
Reconstruct a string from a sequence of (k,d)-mers corresponding to a path in a paired de Bruijn graph.

Given: A sequence of (k, d)-mers (a1|b1), ... , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for all i from 1 to n-1.

Return: A string Text where the i-th k-mer in Text is equal to Suffix(ai|bi) for all i from 1 to n, if such a string exists.

Sample Dataset
4 2
GACC|GCGC
ACCG|CGCC
CCGA|GCCG
CGAG|CCGG
GAGC|CGGA

Sample Output
GACCGAGCGCCGGA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> Construct De Bruijn Graph (BA3D)
      -> Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> Eulerian Cycle (BA3F)
        -> PREV: k-universal string problem (BA3I)
      -> Eulerian Path (BA3G)
      ↓
    * Read-Pair (BA3J)
                                         "
              not every Eulerian path in the paired de Bruijn graph
         constructed from a (k, d)-mer composition spells out a solution of
                 the String Reconstruction from Read-Pairs Problem.
                                         "
      -> HERE: Construct a String Spelled by a Gapped Genome Path (BA3L)
      -> NEXT: Generate All Maximal Non-Branching Paths in a Graph (BA3M)

Info.
  * figure 3.36

            AG    A    CA
       AG   TG ↙ T ↖ CT   CT
       AG   ↙        ↖    CA
     A -> G  ⇛  GC ⇛    C  ->  T
     A    G      GC       C      A
            ↖         ↙
          TG  ↖  T  ↙  CT
          TG      T      CT

  (1) Eulerian Path 1 (incorrect)
    (AG|AG) -> (GC|GC) -> (CT|CT) -> (TG|TG) -> (GC|GC) -> (CA|CT) -> (AG|TG) -> (GC|GC) -> (CT|CA)
    (AG|AG)     A
     (GC|GC)    |  mismatch
      (CT|CT)   T
       (TG|TG)     T
        (GC|GC)    |  mismatch
         (CA|CT)   A
          (AG|TG)
           (GC|GC)
            (CT|CA)
     AGCTGCAGCT
           AGCTGCTGCA    <-- The overlapping part depends on the path.
     AGCTGCAGCTGCTGCA

  (2) Eulerian Path 2 (correct)
    (AG|AG) -> (GC|GC) -> (CA|CT) -> (AG|TG) -> (GC|GC) -> (CT|CT) -> (TG|TG) -> (GC|GC) -> (CT|CA)
    (AG|AG)
     (GC|GC)
      (CA|CT)
       (AG|TG)
        (GC|GC)
         (CT|CT)
          (TG|TG)
           (GC|GC)
            (CT|CA)
     AGCAGCTGCT
        AGCTGCTGCA     <-- The overlapping part depends on the path.
     AGCAGCTGCTGCAGCA
        |  |     |
     AGCTGCAGCTGCTGCA  <-- string in case (1)
     => Not all Eulerian Path in the paired De Bruijn graph
        spell solutions of String Reconstruction from Read-Pairs Problem!

  ╔═════════════════════════════════════════════════════════════════════════════════════╗
  ║ STRINGSPELLEDBYGAPPEDPATTERNS(GappedPatterns, k, d)                                 ║
  ║   FirstPatterns <- the sequence of initial k-mers from GappedPatterns               ║
  ║   SecondPatterns <- the sequence of terminal k-mers from GappedPatterns             ║
  ║   PrefixString <- STRINGSPELLEDBYPATTERNS(FirstPatterns, k)                         ║
  ║   SuffixString <- STRINGSPELLEDBYPATTERNS(SecondPatterns, k)                        ║
  ║   for i = k + d + 1 to |PrefixString|                                               ║
  ║     if the i-th symbol in PrefixString != the (i - k - d)-th symbol in SuffixString ║
  ║       return "there is no string spelled by the gapped patterns"                    ║
  ║   return PrefixString concatenated with the last k + d symbols of SuffixString      ║
  ╚═════════════════════════════════════════════════════════════════════════════════════╝

"""

import time


def str_reconst(path):
    """
    [(a,a)] -> a
    returns a reconstructed string from the result of 'PathGraph_(k,d)(TGext)'
    correct sometimes, not always
    TODO: Find another algorithm. This problem is not as simple as I thought.
    """
    string1 = path[0][0]
    string2 = path[0][1]
    for tup in path[1:]:
        string1 += tup[0][-1]
        string2 += tup[1][-1]
    idx = 0
    for i in range(len(string1)):
        if not string2.find(string1[i:]) == -1:
            idx = i
            break
    return string1[:idx] + string2


def string_spelled_by_patterns(patterns, k):
    """
    ([str],int) -> str
    implementation of 'StringSpelledByPatterns(Patterns, k)'
    """
    result = patterns[0]
    for pattern in patterns[1:]:
      result += pattern[-1]
    return result


def string_spelled_by_gapped_patterns(gpatterns, k, d):
    """
    [(str,str)] -> str
    implementation of 'STRINGSPELLEDBYGAPPEDPATTERNS(GappedPatterns, k, d)'
    TODO: Is that all? This is not the perfect solution I thought of.
    """
    patterns_fst = []
    patterns_snd = []
    for tup in gpatterns:
        patterns_fst.append(tup[0])
        patterns_snd.append(tup[1])
    str_prefix = string_spelled_by_patterns(patterns_fst, k)
    str_suffix = string_spelled_by_patterns(patterns_snd, k)
    for i in range(k+d, len(str_prefix)):
        if not str_prefix[i] == str_suffix[i-k-d]:
            return "There is no string spelled by the gapped patterns"
    return str_prefix + str_suffix[-(k+d):]


def main():
    f = open('/home/wsl/rosalind/data/ba03l.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k, d = list(map(int, lines[0].split()))
    reads = []
    for line in lines[1:]:
        pattern1, pattern2 = [pattern.strip() for pattern in line.split('|')]
        reads.append((pattern1, pattern2))
    f.close()

    start_time = time.time()
    print(str_reconst(reads))
    print(string_spelled_by_gapped_patterns(reads,k,d))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()