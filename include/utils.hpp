/*
utility functions

* read_lines
* split_str
* trim
* read_graph

* print_arr
* print_vec
* print_set
* print_unordered_set
* print_unordered_set_simple
*/

#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>  // find_if
#include <iostream>

// read lines
std::vector<std::string> read_lines(std::string filename) {
    std::ifstream file(filename);
    std::vector<std::string> lines;
    std::string str;
    while (getline(file, str)) {
        lines.push_back(str);
    }
    file.close();
    return lines;
}

// split string
std::vector<std::string> split_str(std::string str, std::string delimiter) {
    std::vector<std::string> splitted;
    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
        splitted.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    splitted.push_back(str);
    return splitted;
}

// trim from start (in place)
// https://stackoverflow.com/questions/216823/how-to-trim-a-stdstring
static inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
        }));
}

static inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
        }).base(), s.end());
}

static inline void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

// read graph
std::unordered_map<std::string, std::vector<std::string>> read_graph(std::vector<std::string> lines) {
    std::unordered_map<std::string, std::vector<std::string>> map;
    for (auto line : lines) {
        std::vector<std::string> parts = split_str(line, " -> ");
        std::string key = parts[0];
        std::vector<std::string> values = split_str(parts[1], ",");
        map[key] = values;
    }
    return map;
}

template <class InputIterator>
void print_arr(InputIterator start, InputIterator end) {
    for (auto itr = start; itr != end; ++itr) {
        std::cout << *itr << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_vec(const T& vec) {
    for (auto iter = vec.begin(); iter != vec.end(); iter++) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_set(const T& st) {
    for (auto iter = st.begin(); iter != st.end(); iter++) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}

// print unordered_map
void print_unordered_map(std::unordered_map<std::string, std::vector<std::string>> map) {
    for (const auto& [key, values] : map) {
        std::cout << key << " -> ";
        std::string val;
        for (auto it = values.begin(); it != values.end(); it++) {
            if (it != values.begin()) {
                val += ",";
            }
            val += *it;
        }
        std::cout << val << std::endl;
    }
}

template <typename T, typename S>
void print_unordered_map_simple(std::unordered_map<T, S> map) {
    for (auto it = map.begin(); it != map.end(); ++it) {
        std::cout << it->first << ": " << it->second << std::endl;
    }
}

#endif