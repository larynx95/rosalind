/*

*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <random>
#include <time.h>

using namespace std;

/************************************************
prototypes
************************************************/

vector<string> read_lines(string filename);
vector<string> split_str(string str, string delimiter);


/************************************************
main
************************************************/

int main() {
    vector<string> lines = read_lines("/home/wsl/rosalind/data/ba03m.txt");

    clock_t tStart = clock();
    printf("Time taken: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

/************************************************
definitions
************************************************/

// helper fx.
vector<string> read_lines(string filename) {
    ifstream file(filename);
    vector<string> lines;
    string str;
    while (getline(file, str))
    {
        lines.push_back(str);
    }
    file.close();
    return lines;
}

// helper fx
vector<string> split_str(string str, string delimiter) {
    vector<string> splitted;
    size_t pos = 0;
    string token;
    while ((pos = str.find(delimiter)) != string::npos) {
        token = str.substr(0, pos);
        splitted.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    splitted.push_back(str);
    return splitted;
}