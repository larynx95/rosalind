"""
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

Return:DeBruijn_k(Text), in the form of an adjacency list.

Sample Dataset
4
AAGATTCTCTAC

Sample Output (1:many)
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC
TTC -> TCT

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * Construct OverapGraph (BA3C)
      -> PREV: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> HERE: Construct De Bruijn Graph (BA3D)
      -> NEXT: Construct De Bruijn Graph with k-mers (BA3E)

Example (Overlap graph):
  TAATGCCATGGGATGTT
  TAA AAT ATG TGC GCC CCA CAT ATG TGG GGG GGA GAT ATG TGT GTT

  Example (De Bruijn graph):
  edges:  TAA AAT ATG TGC GCC CCA CAT ATG TGG GGG GGA GAT ATG TGT GTT
  nodes: TA->AA->AT->TG->GC->CC->CA->AT->TG->GG->GG->GA->AT->TG->GT->TT
  repeat:         *   ^               *   ^   v   v       *   ^

Plan 1.
- AAGATTCTCTAC
  => k-mers    :  AAGA    AGAT    GATT    ATTC    TTCT    TCTC    CTCT    TCTA    CTAC
                  /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \
  => (k-1)-mers: AAG AGA AGA GAT GAT ATT ATT TTC TTC TCT TCT CTC CTC TCT TCT CTA CTA TAC
                 1   2   2   2   2   2   2   2   2   4   4   2   2   4   4   2   2   1
  => (AAG,AGA) (AGA,GAT) (GAT,ATT) (ATT,TTC) (TTC,TCT) (TCT,CTC) (CTC,TCT) (TCT,CTA) (CTA,TAC)
  (AAG,AGA)
      (AGA,GAT)
          (GAT,ATT)
              (ATT,TTC)
                  (TTC,TCT)
                      (TCT,CTC)
                          (CTC,TCT)
                      (TCT,CTA)       <-- TCT again!
                          (CTA,TAC)

═════════════════════════════════════════════════

References:
- How To Pronounce "De Bruijn"?
  https://www.biostars.org/p/7186/
- Python update a key in dict if it doesn't exist
  https://stackoverflow.com/questions/42315072/python-update-a-key-in-dict-if-it-doesnt-exist
- How can I sort a dictionary by key?
  https://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key
- Writing to a File with Python's print() Function
  https://stackabuse.com/writing-to-a-file-with-pythons-print-function/
- Print string to text file
  https://stackoverflow.com/questions/5214578/print-string-to-text-file
- How to redirect 'print' output to a file?
  https://stackoverflow.com/questions/7152762/how-to-redirect-print-output-to-a-file
"""

import time
import sys


def all_kmers(text, k):
    """
    (str,int) -> gen
    returns a generator of all k-mers from a text string
    """
    for i in range(len(text)-k+1):
        yield text[i:i+k]


def de_bruijn(text, k):
    """
    (str,int) -> {str:[str]}
    Implementation of "DeBruijn_k(Text)"
    >>> de_bruijn('AAGATTCTCTAC', 4)
        {'AAG':['AGA'],'AGA':['GAT'],'GAT':['ATT'],'ATT':['TTC'],'TTC':['TCT'],'TCT':['CTC','CTA'],'CTC':['TCT'],'CTA':['TAC']}
    """
    dict = {}
    for kmer in all_kmers(text, k):
        key = kmer[:-1]
        val = kmer[1:]
        if key not in dict:
            dict[key] = [val]
        else:
            dict[key].append(val)
    return dict


def print_de_bruijn(dict):
    sorted_keys = sorted(dict.keys())
    for key in sorted_keys:
        print(key, end=" -> ")
        print(*dict[key], sep =",")


def main():
    f = open('/home/wsl/rosalind/data/ba03d.txt', 'r')
    lines = [line.strip() for line in f.readlines()]
    k = int(lines[0])
    text = lines[1]
    f.close()

    # check time
    start_time = time.time()
    print_de_bruijn(de_bruijn(text, k))  # too long to copy all, create an output file!
    print("--- %s seconds ---" % (time.time() - start_time))

    # print to output file
    start_time = time.time()
    original_stdout = sys.stdout
    with open('/home/wsl/rosalind/data/ba03d_output.txt', 'w') as f:
        sys.stdout = f
        print_de_bruijn(de_bruijn(text, k))
        sys.stdout = original_stdout
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()