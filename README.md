# MPCE
maximal k-partite clique enumaration

## Usage
```
Usage: ./mpce Graph <options>
Options:
-p				print out clique list
	
				note: no print out if not specify -p
	
-o <filename>	filename to store cliques if choose to print out
				print to stdout if not specify -o
-l <value>		least number of vertices in a clique <default = 3>
-u <value>		most number of vertices in a clique <default = graph size>
-v [1~8]	    algorithm version <default = 2>. 4~8 are maximal k-partite clique enumeration algorithms.
					1 - bron kerbosch version 1 (numerical order)
					2 - bron kerbosch version 2 (improved version)
					3 - modified bron kerbosch to find maximum clique
					4 - maximal partite cliques enumration (BK-style k-partite search with pivot). Terminates when the candidate set of a partite is empty and no vertex from that partite is included in the current k-partite clique.
					5 - maximal partite cliques enumration (MMCE add intra-partite edges + BK with pivot/fixp selection)
					6 - maximal partite cliques enumration (MPCE). Terminates when the candidate set from one partite is empty and no vertex from this partite is present in the current clique; also forces the first k selected vertices to come from different partite sets.
					7 - maximal partite cliques enumration (MMCE add intra-partite edges + BK with pivot)
					8 - maximal partite cliques enumration (MMCE variant + BK with pivot)
```
## Example Usage
- Maximal k-partite clique enumration
    - ./mpce test/test_K222.txt -v 6 -klb 1 -p
- The 'make test version_number' command will execute all .txt files located within the test directory using the provided version number. Support for versions 4 to 8.  

## Input file format requirement when enumerate maximal k-partite cliques (version 4- 8)
- Comment line
	- The first row should begin with the character `c` followed by a space.  
	- This line is treated as a **comment** or **explanation** of the file's contents.
- Graph matrics line
	- The second line must begin with "p edges", followed by three integers: vertex_number edge_number partite_set_number
- Edge lines
	- Each row represents an edge, using the format "e node1 node2"
- Node lines
	- Each row represents a node, using the format "n node partite"

## Implementation Notes
- `v5` does use pivoting. Its `fixp` selection is the classic Bron-Kerbosch/Tomita-style rule: choose a pivot from `X union P` that maximizes the number of neighbors in `P` (equivalently, minimizes the number of non-neighbors in `P`).
- The shared Tomita pivot helpers were corrected to search pivot candidates in `X union P`, while still scoring them only by their compatibility with vertices in `P`.
- For the k-partite routines, same-partite checks now use a bitset-backed `same_category(...)` lookup, similar to the bitset-backed `edge_exists(...)` test.
- A cleanup bug that could abort after printing without `-o` was fixed by avoiding `fclose(stdout)`.
