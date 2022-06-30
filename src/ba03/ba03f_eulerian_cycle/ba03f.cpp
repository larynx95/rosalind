/*
Rosalind: BA3F (difficulty: 4/5)
Find an Eulerian Cycle in a Graph

A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle,
and we say that a graph containing such a cycle is Eulerian.
The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

╔══════════════════════════════════════════════════════════════════════╗
║ EULERIANCYCLE(Graph)                                                 ║
║     form a cycle Cycle by randomly walking in Graph                  ║
║         (don't visit the same edge twice!)                           ║
║     while there are unexplored edges in Graph                        ║
║         select a node newStart in Cycle with still unexplored edges  ║
║         form Cycle' by traversing Cycle (starting at newStart)       ║
║           and then randomly walking                                  ║
║         Cycle <- Cycle'                                              ║
║     return Cycle                                                     ║
╚══════════════════════════════════════════════════════════════════════╝

╔═════════════════════════════════════════════════════════════╗
║ ALLEULERIANCYCLES(Graph)                                    ║
║     AllGraph <- the set consisting of a single graph Graph  ║
║     while there is a non-simple graph G in AllGraphs        ║
║         v <- a node with indegree larger than 1 in G        ║
║         for each incoming edge (u,v) into v                 ║
║             for each outgoing edge (v,w) from v             ║
║                 NewGraph <- (u,v,w)-bypass graph of G       ║
║                 if NewGraph is connected                    ║
║                     add NewGraph to ALlGraphs               ║
║         remove G from AllGraphs                             ║
║     for each graph G in AllGraphs                           ║
║         output the (single) Eulerian cycle in G             ║
╚═════════════════════════════════════════════════════════════╝

Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An "Eulerian directed graph", in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3
1 -> 0
2 -> 1,6
3 -> 2
4 -> 2
5 -> 4
6 -> 5,8
7 -> 9
8 -> 7
9 -> 6

Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

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
      -> PREV: Construct De Bruijn Graph with k-mers (BA3E)
      ↓
    * Eulerian Cycle, Path
      -> HERE: Eulerian Cycle (BA3F)
      -> NEXT: Eulerian Path (BA3G)

Info.
1. Hierholzer's algorithm
  * steps:
    ; just walk through the graph, and get a circle1
    ; if the circle1 is complete whole large cycle without remaining nodes in graph, that's Eulerian cycle
    ; if there're remaining unvisited node in graph yet, just do another walking, get another cycle2
      but in this time, start node should be a node in previous cycle, and it has one or more out-degrees
    ; patch cycle2 to cycle1, repeat previous steps
              4 < 5             4 < 5             4 * 5             4 * 5
              v   ^             v   ^             *   *             *   *
          1 < 2 > 6 > 8     1 *[2]> 6 > 8     1 * 2 * 6 > 8     1 * 2 *[6]* 8
          v   ^   ^   v     *   *   ^   v     *   *   ^   v     *   *   *   *
          0 > 3   9 < 7     0 * 3   9 < 7     0 * 3   9 < 7     0 * 3   9 * 7
                            10321             10321             103265421
                                                 26542              68796
                                              103265421         1032687965421

═════════════════════════════════════════════════

References:
-
*/

#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>
#include "/home/wsl/rosalind/include/utils.hpp"  // includes utility functions

using namespace std;
// typedef std::unordered_map<std::string, std::vector<std::string>> umap;   // typedef

/************************************************
prototypes
************************************************/

vector<string> walk(unordered_map<string, vector<string>>& graph, string node);
bool is_depleted(unordered_map<string, vector<string>>& graph);
vector<string> eulerian_cycle(unordered_map<string, vector<string>>& graph);
void print_eulerian_cycle(std::vector<std::string> const& cycle);

/************************************************
main
************************************************/

int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03f.txt");

    clock_t tStart = clock();
    unordered_map<string, vector<string>> graph = read_graph(lines);
    vector<string> cycle = eulerian_cycle(graph);
    print_eulerian_cycle(cycle);
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

/************************************************
definitions
************************************************/

// simple walk
// graph should be passed by reference, because it must be modified.
// C++ support both passed-by-value and passed-by-reference.
vector<string> walk(unordered_map<string, vector<string>>& graph, string source) {
    srand(time(NULL));
    vector<string> cycle = { source };
    while (graph.at(source).size() != 0) {
        // random select target node
        vector<string> values = graph[source];
        int random_idx = rand() % values.size();
        string target = values[random_idx];
        // remove the target node from graph[node]
        values.erase(values.begin() + random_idx);
        graph.at(source) = values;
        // add target node to cycle
        cycle.push_back(target);
        // update source node
        source = target;
    }
    return cycle;
}

// check if all values are empty
bool is_depleted(unordered_map<string, vector<string>>& graph) {
    for (auto& [key, values] : graph) {
        if (!(values.empty())) return false;
    }
    return true;
}

// print Eulerian cycle
void print_eulerian_cycle(std::vector<std::string> const& cycle) {
    for (auto it = cycle.begin(); it != cycle.end(); it++) {
        if (it == cycle.end() - 1) {
            std::cout << *it << std::endl;
            break;
        } else {
            std::cout << *it << "->";
        }
    }
}

// BA03F: Eulerian cycle
vector<string> eulerian_cycle(unordered_map<string, vector<string>>& graph) {
    // get all keys from a map (unordered_map)
    vector<string> nodes;
    for (auto& pair : graph) {
        nodes.push_back(pair.first);
    }
    // select a start node from the vector of nodes
    srand(time(NULL));
    string start = nodes[rand() % nodes.size()];
    // get a cycle by simple walking (main cycle which many new cycles will be merged into)
    vector<string> cycle = walk(graph, start);
    // iteration
    while (!is_depleted(graph)) {
        // select a node to get another cycle
        for (auto node : cycle) {
            if (graph[node].size() != 0) {
                start = node;
                break;
            }
        }
        // get the index of start node in the previous (main) cycle
        int idx = 0;
        for (size_t i = 0; i < cycle.size(); i++) {
            if (cycle[i] == start) {
                idx = i;
                break;
            }
        }
        // get new cycle
        vector<string> new_cycle = walk(graph, start);
        // merging new cyle to previous cycle
        new_cycle.pop_back();
        cycle.insert(cycle.begin() + idx, new_cycle.begin(), new_cycle.end());
    }
    return cycle;
}
