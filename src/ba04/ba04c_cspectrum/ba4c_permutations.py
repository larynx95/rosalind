'''
Permutations

1. Treminology:
permutation (order): how to arrange n objects in n positions
combination (inorder)
cartesian
product

2. Details
1) permutation (order): how to arrange n objects in n positions
- the number of permutations and combinations
  nCr = n! / r! (n - r)!
  nPr = n! / (n - r)!
- Pseudocode for Simple Recursive Algorithm
    Data:
    GeneratedPermutations: generated prmutations
    CurrentPermutation: current permutation
    ElementsToPermute: elements to permute
    if ElementToPermute not empty then
        for element in ElementsTopermute do
        NextPermutation <- Currentpermutation + element;
        RemainingElements <- remove element from ElementsToPermute
        Permutation(GeenratedPermutations, NextPermutation, RemainingElements);
        end
    else
        add CurrentPermutation to GeneratedPermutations;
    end
- Pseudocode for Heap's algorithm
  This algorithm is based on swapping elements to generate the permutations.
  The principle of Heap’s algorithm is decrease and conquer.
  https://www.geeksforgeeks.org/decrease-and-conquer/

2) combination (inorder)
- the number of permutations and combinations
  nCr = n! / r! (n - r)!
  nPr = n! / (n - r)!

3) cartesian product

References:
- Generate All Permutations of an Array
  https://www.baeldung.com/cs/array-generate-all-permutations
'''

def all_combs_rec(text):
    '''
    >>> all_combs_rec('NQE')
        ['', 'N', 'Q', 'QN', 'E', 'EN', 'EQ', 'EQN']
    '''
    if len(text) == 0:
        return ['']
    cs = []
    for c in all_combs_rec(text[1:]):
        cs += [c, c+text[0]]
    return cs


def all_combs_gen(text):
    '''
    all combinations, generator version
    >>> list(all_combs('NQE'))
        ['', 'N', 'Q', 'QN', 'E', 'EN', 'EQ', 'EQN']
    >>> print(*all_combs_gen('LQNE'))
        L Q QL N NL NQ NQL E EL EQ EQL EN ENL ENQ ENQL
    >>> print(*all_combs_gen('NQEL'))
        N Q QN E EN EQ EQN L LN LQ LQN LE LEN LEQ LEQN
    '''
    if len(text) == 0:
        yield ''
    else:
        for c in all_combs_gen(text[1:]):
            yield c
            yield c+text[0]


def all_perms(text):
    '''
    >>> list(all_perms('NQE'))
        ['NQE', 'QNE', 'QEN', 'NEQ', 'ENQ', 'EQN']
    '''
    if len(text) <=1:
        yield text
    else:
        for perm in all_perms(text[1:]):
            for i in range(len(text)):
                yield perm[:i] + text[0:1] + perm[i:]