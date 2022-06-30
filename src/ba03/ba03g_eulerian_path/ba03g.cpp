/*
Rosalind: BA3G
Find an Eulerian Path in a Graph

In "Find an Eulerian Cycle in a Graph", we defined an Eulerian cycle.
A path that traverses each edge of a graph exactly once
(but does not necessarily return to its starting node is called an Eulerian
path.

Eulerian Path Problem
Find an Eulerian path in a graph.

Given: A directed graph that contains an Eulerian path,
       where the graph is given in the form of an adjacency list.

Return: An Eulerian path in this graph.

Sample Dataset
0 -> 2
1 -> 3
2 -> 1
3 -> 0,4
6 -> 3,7
7 -> 8
8 -> 9
9 -> 6

Sample Output
6->7->8->9->6->3->0->2->1->3->4

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
      -> PREV: Eulerian Cycle (BA3F)
        -> NEXT: k-universal string problem (BA3I)
      -> HERE: Eulerian Path (BA3G)

Info.
  * writing down graph to figure out what this exercise is:
    {0:[2], 1:[3], 2:[1], 3:[0, 4], 6:[3,7], 7:[8], 8:[9], 9:[6]}
                  4
                  ↑
          9 → 6 → 3 → 0       Not circular!!!
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
  * problems:
    - This graph is not circular graph.
      Can any node be a starting point (source) or an ending point (target)?
    - How can I find the end node?
      a. comparing keys and values in dictionary (De Bruijn graph)
         dict.keys:   0123  6789
         dict.values: 012334 789   --> 4 is the end!
      b. difference bw in-degree and out-degree
         In-degree must be larger than out-degree in the last node.
         0:0, 1:0, 2:0, 3:0, 4:1, 6:-1, 7:0, 8:0, 9:0  --> The node 4 is the
last!

Plan 1.
  * re-using "EulerianCycle" function
    - Adding an extra edge from w (end node) to v (start node)
      transforms the Eulerian path into an Eulerian cycle.
    - convert Eulerian path to Eulerian cycle, then re-convert cycle to path

  * out-degree +1, in-degree -1
                  4
                  ↑
          9 → 6 → 3 → 0
          ↑   ↓   ↑   ↓
          8 ← 7   1 ← 2
    unique nodes: 0    1    2    3    4    6    7    8    9
    --------------------------------------------------------
      0:[2]       +1       -1
      1:[3]           +1        -1
      2:[1]           -1   +1
      3:[0,4]     -1            +2   -1
      6:[3,7]                   -1        +2   -1
      7:[8]                                    +1   -1
      8:[9]                                         +1   -1
      9:[6]                               -1             +1
    --------------------------------------------------------
                  0    0    0    0   -1   +1    0    0    0   --> start:6 end
(sink):4
  * add edge (end, start) to graph, get Eulerian cycle: 4 -> 6
  * trim the last repetitive elem in cycle, then rotate cycle to get Eulerian
path

═════════════════════════════════════════════════

References:
- tuple
  https://en.cppreference.com/w/cpp/utility/tuple
- vector::at vs. vector::operator[]
  https://stackoverflow.com/questions/9376049/vectorat-vs-vectoroperator
*/

#include <time.h>

#include <limits>
#include <random>
#include <unordered_set>
#include <tuple>

#include "/home/wsl/rosalind/include/utils.hpp"

using namespace std;

/************************************************
prototypes
************************************************/

vector<string> walk(unordered_map<string, vector<string>>& graph, string node);
bool is_depleted(unordered_map<string, vector<string>>& graph);
vector<string> eulerian_cycle(unordered_map<string, vector<string>>& graph);
unordered_map<string, vector<int>> get_degrees(
    unordered_map<string, vector<string>> const& graph);
tuple<string, string> find_start_sink(unordered_map<string, vector<string>> const& graph);
vector<string> eulerian_path(unordered_map<string, vector<string>>& graph);
void print_degrees(std::unordered_map<std::string, std::vector<int>> const& degrees);
void print_vector_with_arrow(std::vector<std::string> const& cycle);

/************************************************
main
************************************************/

int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03g.txt");
    unordered_map<string, vector<string>> graph = read_graph(lines);

    clock_t tStart = clock();
    unordered_map<string, vector<int>> dgs = get_degrees(graph);
    vector<string> path = eulerian_path(graph);
    print_vector_with_arrow(path);

    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

/************************************************
definitions
************************************************/

// simple walk
vector<string> walk(unordered_map<string, vector<string>>& graph,
    string source) {
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
    // get a cycle by simple walking (main cycle which many new cycles will be
    // merged into)
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

// get degrees of a graph
unordered_map<string, vector<int>> get_degrees(
    unordered_map<string, vector<string>> const& graph) {
    unordered_map<string, vector<int>> degrees;   // unordered_map which will be returned
    for (const auto& [key, values] : graph) {
        // out-degrees
        if (degrees.find(key) == degrees.end()) {
            degrees[key] = { 0, (int)values.size() };
        } else {
            // degrees[key][1] += 1;
            degrees[key][1] += values.size();
        }
        // in-degrees
        for (auto val : values) {
            if (degrees.find(val) == degrees.end()) {
                degrees[val] = { 1,0 };
            } else {
                degrees[val][0] += 1;
            }
        }
    }
    return degrees;
}

// get degrees of a graph
tuple<string, string> find_start_sink(
    unordered_map<string, vector<string>> const& graph) {
    unordered_map<string, vector<int>> degrees;   // unordered_map which will be returned
    for (const auto& [key, values] : graph) {
        // out-degrees
        if (degrees.find(key) == degrees.end()) {
            degrees[key] = { 0, (int)values.size() };
        } else {
            // degrees[key][1] += 1;
            degrees[key][1] += values.size();
        }
        // in-degrees
        for (auto val : values) {
            if (degrees.find(val) == degrees.end()) {
                degrees[val] = { 1,0 };
            } else {
                degrees[val][0] += 1;
            }
        }
    }
    string start, sink;
    for (auto const& [key, values] : degrees) {
        if (values[1] - values[0] > 0) {  // if out-degrees > in-degrees
            start = key;
        } else if (values[1] - values[0] < 0) {  // if out-degrees < in-degrees
            sink = key;
        }
    }
    return make_tuple(start, sink);
}

// print degrees
void print_degrees(std::unordered_map<std::string, std::vector<int>> const& degrees) {
    for (auto const& [key, values] : degrees) {
        std::cout << key << ":";
        for (size_t i = 0; i < values.size(); i++) {
            if (i == 1) {
                std::cout << values[i] << std::endl;
                break;
            } else {
                std::cout << values[i] << ",";
            }
        }
    }
}

// print Eulerian cycle
void print_vector_with_arrow(std::vector<std::string> const& cycle) {
    for (auto it = cycle.begin(); it != cycle.end(); it++) {
        if (it == cycle.end() - 1) {
            std::cout << *it << std::endl;
            break;
        } else {
            std::cout << *it << "->";
        }
    }
}

// BA03G: Eulerian path
vector<string> eulerian_path(unordered_map<string, vector<string>>& graph) {
    // find start, sink nodes
    tuple<string, string> tup = find_start_sink(graph);
    string start = get<0>(tup);
    string sink = get<1>(tup);
    // connect start and sink nodes: path -> cycle
    graph[sink].push_back(start);  // TODO: graph.at(sink).push_back(start);  <-- This is wrong, why?
    vector<string> cycle = eulerian_cycle(graph);
    // find index of break point (index of sink node)
    auto it = find(cycle.cbegin(), cycle.cend(), sink);
    int idx = it - cycle.cbegin();
    // disconnect, transposition
    vector<string> path = { cycle.begin() + idx + 1, cycle.end() - 1 };
    vector<string> rest = { cycle.begin(), cycle.begin() + idx + 1 };
    path.insert(path.end(), rest.begin(), rest.end());
    return path;
}
