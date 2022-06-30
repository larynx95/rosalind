/*
Rosalind: BA3L
Construct a String Spelled by a Gapped Genome Path

Gapped Genome Path String Problem
Reconstruct a string from a sequence of (k,d)-mers corresponding to a path in a paired de Bruijn graph.

Given:
A sequence of (k, d)-mers (a1|b1), ... , (an|bn)
such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for all i from 1 to n-1.

Return:
A string Text where the i-th k-mer in Text is equal to
Suffix(ai|bi) for all i from 1 to n, if such a string exists.

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

  * sample dataset
    GACC--GCGC
     ACCG--CGCC
      CCGA--GCCG
       CGAG--CCGG
        GAGC--CGGA
    GACCGAGC
          GCGCCGGA
    GACCGAGCGCCGGA

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
*/

// BA03L: string spelled by patterns
function str_spelled_by_patterns(patterns, k) {
    /**
     * ([str],int) -> str
     * returns a string from FirstPatterns
     */
    let result = patterns[0];
    for (const pat of patterns.slice(1)) {
        result += pat[pat.length - 1];
    }
    return result;
}

// BA03L: string spelled by gapped patterns
function str_spelled_by_gapped_patterns(gpatterns, k, d) {
    /**
     * ([str],int,int) -> str
     * returns a string from gapped patterns
    */
    let pat_fst = [];
    let pat_snd = [];
    for (const gpat of gpatterns) {
        const [fst, snd] = gpat.split('|');
        pat_fst.push(fst);
        pat_snd.push(snd);
    }
    str_prefix = str_spelled_by_patterns(pat_fst, k);
    str_suffix = str_spelled_by_patterns(pat_snd, k);
    for (let i = k + d; i < str_prefix.length; i++) {
        if (str_prefix[i] !== str_suffix[i - k - d]) return "No stirng!";
    }
    return str_prefix + str_suffix.slice(str_suffix.length - (k + d));
}

// main
function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba03l.txt').toString().split("\n");
    const [k, d] = lines[0].split(/\s/).map(Number);  // Don't forget parseInt()!
    const gpatterns = lines.slice(1);

    const startTime = performance.now()
    const gstr = str_spelled_by_gapped_patterns(gpatterns, k, d);
    console.log(gstr);
    console.log(`${performance.now() - startTime} milliseconds`)
}

// execute main function
main()