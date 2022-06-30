# Rosalind

## Chapter 03. How Do We Assemble Genomes? (graph algorithm)
- [BA3A. Generate the k-mer Composition of a String1](https://rosalind.info/problems/ba3a/)
- [BA3B. Reconstruct a String from its Genome Path](https://rosalind.info/problems/ba3b/)
- [BA3C. Construct the Overlap Graph of a Collection of k-mers](https://rosalind.info/problems/ba3c/)
- [BA3D. Construct the De Bruijn Graph of a String](https://rosalind.info/problems/ba3d/)
- [BA3E. Construct the De Bruijn Graph of a Collection of k-mers](https://rosalind.info/problems/ba3e/)
- **[BA3F. Find an Eulerian Cycle in a Graph](https://rosalind.info/problems/ba3f/)**
- [BA3G. Find an Eulerian Path in a Graph](https://rosalind.info/problems/ba3g/)
- [BA3H. Reconstruct a String from its k-mer Composition](https://rosalind.info/problems/ba3h/)
- **[BA3I. Find a k-Universal Circular String](https://rosalind.info/problems/ba3i/)**
- [BA3J. Reconstruct a String from its Paired Composition](https://rosalind.info/problems/ba3j/)
- **[BA3K. Generate Contigs from a Collection of Reads](https://rosalind.info/problems/ba3k/)**
- [BA3L. Construct a String Spelled by a Gapped Genome Path](https://rosalind.info/problems/ba3l/)
- **[BA3M. Generate All Maximal Non-Branching Paths in a Graph](https://rosalind.info/problems/ba3m/)**

---
## ☘ Info
- Eulerian cycle
  - Hierholzer's algorithm
  - Fleury's Algorithm
  - [Eunice(?)'s algorithm](https://math.stackexchange.com/questions/3965493/determine-if-edge-is-a-bridge-in-a-graph)
  - ratating cycle
  - recursive, closure
  - Tree
  - merging small cycles to large one

---
## ☘ Todo
- (BA03D, BA03D) Is duplication allowed in De Bruijn graph? Which is correct?
    ```go
    // BA03D: De Bruijn graph from a text
    func DeBruijnFromText(text string, k int) map[string][]string {
    	graph := make(map[string][]string)
    	for i := 0; i < len(text)-k+1; i++ {
    		prefix := text[i : i+k-1]
    		suffix := text[i+1 : i+k]
    		if values, ok := graph[prefix]; ok {
    			values = append(values, suffix) // TODO:  can be duplicated
    			graph[prefix] = values
    		} else {
    			graph[prefix] = []string{suffix}
    		}
    	}
    	return graph
    }

    // BA03D: De Bruijn graph from a text, set version
    func DeBruijnFromText_set(text string, k int) map[string](map[string]int) {
    	graph := make(map[string]map[string]int)
    	for i := 0; i < len(text)-k+1; i++ {
    		prefix := text[i : i+k-1]
    		suffix := text[i+1 : i+k]
    		if _, ok := graph[prefix]; ok {
    			graph[prefix][suffix] = 1
    		} else {
    			graph[prefix] = map[string]int{suffix: 1}
    		}
    	}
    	return graph
    }
    /*
    text: "AAGATTCTCTACAAGATTCTCTAC" ("AAGATTCTCTAC" + "AAGATTCTCTAC")
    k: 4
    the first function                the second function
    -----------------------------------------------------
    AAG -> AGA,AGA                    TTC -> TCT
    ATT -> TTC,TTC                    CTA -> TAC
    TTC -> TCT,TCT                    TAC -> ACA
    TAC -> ACA                        ATT -> TTC
    CAA -> AAG                        TCT -> CTC,CTA
    AGA -> GAT,GAT                    CTC -> TCT
    GAT -> ATT,ATT                    ACA -> CAA
    TCT -> CTC,CTA,CTC,CTA            CAA -> AAG
    CTC -> TCT,TCT                    AAG -> AGA
    CTA -> TAC,TAC                    AGA -> GAT
    ACA -> CAA                        GAT -> ATT
    */
    ```

