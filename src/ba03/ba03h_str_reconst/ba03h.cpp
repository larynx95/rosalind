/*
Rosalind: BA3H
Reconstruct a String from its k-mer Composition

String Reconstruction Problem
Reconstruct a string from its k-mer composition.

Given: An integer k followed by a list of k-mers Patterns.

Return: A string Text with k-mer composition equal to Patterns.
(If multiple answers exist, you may return any one.)

Sample Dataset
4
CTTA
ACCA
TACC
GGCT
GCTT
TTAC

Sample Output
GGCTTACCA

═════════════════════════════════════════════════

    [ Where am I? ]

    * How do we assemble genome?
      ↓
    * String reconstruction problem
      -> Generate the k-mer Composition of a String (BA3A)
      -> PREV: String Reconstruction Problem with k-mer composition(BA3B)
      ↓
    * PREV: Construct OverapGraph (BA3C)
      -> HERE: Reconstruct a String from its k-mer Composition (BA3H)
      ↓
    * De Bruijn Graph
      -> NEXT: Construct De Bruijn Graph (BA3D)

Plan 1.
  * multiple strings -> overlap graph -> Hamiltonian path
  * split k-mer to two (k-1)-mers (prefix and suffix)
    CTTA -> CTT,TTA
    ACCA -> ACC,CCA
    TACC -> TAC,ACC
    GGCT -> GGC,GCT
    GCTT -> GCT,CTT
    TTAC -> TTA,TAC
    prefixes: [CCT,ACC,TAC,GGC,GCT,TTA]
    suffixes: [TTA,CCA,ACC,GCT,CTT,TAC]

Plan 2.
  * using Eulerian path (BA3G) algorithm
  * multiple strings -> De Bruijn graph -> Eulerian path

═════════════════════════════════════════════════

References:
-
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

unordered_map<string, vector<string>> debruijn_from_kmers(vector<string> const& patterns);
vector<string> walk(unordered_map<string, vector<string>>& graph, string node);
bool is_depleted(unordered_map<string, vector<string>>& graph);
vector<string> eulerian_cycle(unordered_map<string, vector<string>>& graph);
unordered_map<string, vector<int>> get_degrees(
    unordered_map<string, vector<string>> const& graph);
tuple<string, string> find_start_sink(unordered_map<string, vector<string>> const& graph);
vector<string> eulerian_path(unordered_map<string, vector<string>>& graph);
void print_degrees(std::unordered_map<std::string, std::vector<int>> const& degrees);

/************************************************
main
************************************************/

int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03h.txt");
    int k = stoi(lines[0]);
    vector<string> patterns{ lines.begin() + 1, lines.end() };

    clock_t tStart = clock();
    unordered_map<string, vector<string>> graph = debruijn_from_kmers(patterns);
    vector<string> path = eulerian_path(graph);

    string answer = path[0];
    for (size_t i = 1; i < path.size(); i++) {
        answer += path[i][k - 2];
    }
    cout << answer << endl;

    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

/************************************************
definitions
************************************************/

// BA03E: De Bruijn graph from k-mer patterns
unordered_map<string, vector<string>> debruijn_from_kmers(vector<string> const& patterns) {
    unordered_map<string, vector<string>> map;
    for (auto pattern : patterns) {
        string prefix = pattern.substr(0, pattern.length() - 1);
        string suffix = pattern.substr(1);
        if (map.find(prefix) == map.end()) {  // if key is not in map
            map.insert({ prefix, {suffix} });
        } else { // if key is in map
            map.at(prefix).push_back(suffix);
        }
    }
    return map;
}

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
    // disconnect, transposition, TODO: Think about 'iterator out of index' problem!
    if (cycle[cycle.size() - 1] == sink) {
        return { cycle.begin() + 1, cycle.end() };
    }
    vector<string> path = { cycle.begin() + idx + 1, cycle.end() - 1 };
    vector<string> rest = { cycle.begin(), cycle.begin() + idx + 1 };
    path.insert(path.end(), rest.begin(), rest.end());

    return path;
}
