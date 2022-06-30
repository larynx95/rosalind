"""
get all paths from a graph

Problems:
  * representations of graph without weights
    graph :: {a:[a]}
    graph :: {(a:a)}
  * representations of graph with weights
    graph :: {a:[(a,int)]}
    graph :: {(a,a):int}

Example:
  * input:
      { a:[b],d:[e],k:[i],q:[r],c:[f],g:[h],j:[i],b:[c,e],e:[f,j],f:[g],i:[f],s:[],t:[],m:[n],n:[o,v],o:[p],p:[m],v:[w],w:[x],x:[n] }

      a -> b -> c                p -> m       s      t
           ↓    ↓                ↑    ↓
      d -> e -> f -> g -> h      o <- n -> v
           ↓    ↑                     ↑    ↓
           j -> i <- k                x <- w
                ↓
                l                q -> r

  * all path:
      abcfgh  abefgh  abejifgh  abejiil  defgh  dejifgh  dejil  kifgh  kil
      pmnop   vwxnv   pmnvwxnop
      qr      s       t

  * nodes
    - start nodes    : a, d, k, q
    - isolated nodes : s, t
    - sink nodes     : h, l, r
    - branching nodes: b, e, f, i
    - cycle          : m, n, o, p, v, w, x

  * paths
    - two part: non-cycle, cycle
    (1) non-cycle
      - find nodes without in-degree edges: a, d, k, q, s, t
      - branching path at a node with multiple out-degree edges: b, e, f, i
        branching             paths
        -----------------------------
        ab ┌ cf - gh          abcfgh
           └ e ┌ fgh          abefgh
               └ ji ┌ fgh     abejifgh
                    └ l       abejil
        de ┌ fgh              defgh
           └ ji ┌ fgh         dejifgh
                └ l           dejil
        ki ┌ fgh              kifgh
           └ l                kil
        qr                    qr
        s                     s
        t                     t
      - ideas
        The same thing happens again and again.
        - meet a node with multiple out-degrees
        - then do branching
        But depth of branches varies.
        So I can't solve this exercise with simple loop.
    (2) cycle
      - the rest of nodes

  * non-cyclic directed graph
      s      t

      q -> r                              a b c d e f g h i j k l
                                        a - 1 - - - - - - - - - -
      a -> b -> c                       b - - 1 - 1 - - - - - - -
           ↓    ↓                       c - - - - - 1 - - - - - -
      d -> e -> f -> g -> h             d - - - - 1 - - - - - - -
           ↓    ↑                       e - - - - - 1 1 - - - - -
           j -> i <- k                  f - - - - - - 1 - - - - -
                ↓                       g - - - - - - - 1 - - - -
                l                       h - - - - - - - - - - - -
                                        i - - - - - 1 - - - - - 1
                                        j - - - - - - - - 1 - - -
                                        k - - - - - - - - 1 - - -
                                        l - - - - - - - - - - - -

═════════════════════════════════════════════════

References:
- Deep copy of a dict in python
  https://stackoverflow.com/questions/5105517/deep-copy-of-a-dict-in-python
  copied = copy.deepcopy(original_dictionary)
- Directed acyclic graph
  https://en.wikipedia.org/wiki/Directed_acyclic_graph
- What does the "yield" keyword do?
  https://stackoverflow.com/questions/231767/what-does-the-yield-keyword-do
- Find all paths in DAG from starting node
  https://stackoverflow.com/questions/54441582/find-all-paths-in-dag-from-starting-node

"""

import time
import copy


def main():
    start_time = time.time()
    #graph = {'a':['b'],'d':['e'],'k':['i'],'q':['r'],'c':['f'],'g':['h'],'j':['i'],'b':['c','e'],'e':['f','j'],'f':['g'],'i':['f'],'s':[],'t':[],'m':['n'],'n':['o'],'o':['p'],'p':['m'],'v':['w'],'w':['x'],'x':['n'] }
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
