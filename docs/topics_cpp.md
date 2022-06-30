# Topics - C++

- [Topics - C++](#topics---c)
  - [Numeric data](#numeric-data)
    - [Infinity](#infinity)
    - [Random, Weighted Random](#random-weighted-random)
  - [String](#string)
    - [Trim](#trim)
  - [Collections: Arrays](#collections-arrays)
  - [Collections: Vectors](#collections-vectors)
    - [insert, remove element](#insert-remove-element)
    - [check if empty](#check-if-empty)
    - [slice, subvector](#slice-subvector)
    - [concatenation of two vectors](#concatenation-of-two-vectors)
  - [Collections: Tuple (immutable)](#collections-tuple-immutable)
  - [Collections: `unordered_map`](#collections-unordered_map)
    - [get value by key](#get-value-by-key)
    - [insert, delete element to `unordered_map`](#insert-delete-element-to-unordered_map)
    - [iterating over `unordered_map`](#iterating-over-unordered_map)
    - [check if key exists in map](#check-if-key-exists-in-map)
  - [Collections: Unordered_set](#collections-unordered_set)
  - [IO](#io)
    - [Read text file](#read-text-file)
    - [file access](#file-access)
  - [Project Structure](#project-structure)
    - [Simple mono-file project without project](#simple-mono-file-project-without-project)
    - [seperation `*.h` and `*.cpp` file, why?](#seperation-h-and-cpp-file-why)
    - [Multi-file project](#multi-file-project)
      - [> Project by Visual Studio](#-project-by-visual-studio)
      - [> Project by CMake](#-project-by-cmake)
  - [External Libraries](#external-libraries)
    - [Boost library](#boost-library)
    - [vcpkg](#vcpkg)

## Numeric data

### Infinity

- positive/negative infinity (double): `std::numeric_limits<double>::infinity()`
- positive/negative infinity (int): `std::numeric_limits<int>::max()`

### Random, Weighted Random

- random integer
    ```cpp
    srand(time(NULL));
    for (auto _ : "abc"){
        cout << rand()%6 << endl;
    }
    ```
- weighted random
    ```cpp
    #include <iostream>
    #include <random>
    #include <time.h>
    #include <random>
    #include <vector>

    using namespace std;

    int main() {
        std::vector<double> weights{50,56,4};
        std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
        std::mt19937 gen;
        gen.seed(time(0));//if you want different results from different runs
        int N = 10;
        std::vector<int> samples(N);
        for(auto & i: samples){
            i = dist(gen);
            cout << i << endl;
        }
        return 0;
    }
    ```

---
## String

### Trim

- How to trim a std::string?
- [How to trim a std::string?](https://stackoverflow.com/questions/216823/how-to-trim-a-stdstring)
- [Trim whitespace from a String](https://stackoverflow.com/questions/25829143/trim-whitespace-from-a-string)

    ```cpp
    #include <iostream>
    #include <algorithm>  // find_if
    using namespace std;

    // trim from start (in place)
    static inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    static inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    // trim from both ends (in place)
    static inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }

    int main(void) {
    	string s = "  tes t ";
    	trim(s);
    	cout << s << "!" << endl;
    	return 0;
    }
    ```

---
## Collections: Arrays

- return array of strings
    ```cpp
    string* getNames() {
        string* names = new string[3];  // memory allocation, array declaration needs size
        names[0] = "Simon";
        names[1] = "Peter";
        names[2] = "Dave";

        return names;
    }

    delete[] names;  // free memory
    ```

---
## Collections: Vectors

### insert, remove element
- insert an element to vector
- remove an element of vector by index or value
- get an element of vector by index, and remove it from the vector
    ```cpp
    vector<string> vec = {"aa", "bb", "cc"};
    cout << vec.at(1) << endl;   // at(index) in vector, at(key) in map
    vec.erase(vec.begin() + 1);
    print_vec(vec);
    ```
- [inserting a vector in a certain position in another vector](https://stackoverflow.com/questions/60138069/inserting-a-vector-in-a-certain-position-in-another-vector)
    ```cpp
    vecta.insert(pos, vectb.begin(), vectb.end());
    ```
- meaning of `vec.end()`: [C++ Vector Library - end() Function], returns an iterator which points to past-the-end element in the vector container
    ```cpp
    // BA03F: Eulerain cycle
    vector<string> a = {"0","1","2","3","0"};
    vector<string> b = {"1","7","8","1"};
    cout << *(b.end()-1) << endl; // "1"
    b.erase(b.end());
    a.insert(a.begin()+1, b.begin(), b.end());
    for (auto elem: a) {
        cout << elem << " ";
    }
    cout << endl;
    ```

### check if empty
- check if a vector is empty

    ```cpp
    if (vec.empty()) { ...}
    ```

### slice, subvector
- get subvector: [Best way to extract a subvector from a vector?](https://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector)
    ```cpp
    vector<string> vec = {"aa", "bb", "cc", "dd", "ee", "ff"};
    vector<string> subvec1 = {vec.begin()+1, vec.end()};
    vector<string> subvec2 = {vec.begin()+3, vec.begin()+4};
    print_vec(subvec1);
    print_vec(subvec2);
    ```

### concatenation of two vectors
- [Concatenating two std::vectors](https://stackoverflow.com/questions/201718/concatenating-two-stdvectors)
    ```cpp
    vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
    ```
- reposition of subvectors: abcdefg -> defgabc

---
## Collections: Tuple (immutable)

- [Setting the Values of Tuple types in C++](https://stackoverflow.com/questions/23277099/setting-the-values-of-tuple-types-in-c)
    ```cpp
    std::tuple<int, float> t;
    std::get<0>(t) = 4;             // set int to 4
    std::get<1>(t) = 3.45;          // set float to 3.45
    t = std::make_tuple(4, 3.45);   // set both values
    ```

---
## Collections: `unordered_map`

### get value by key

- get value of a key in `unordered_map`: `map.at(key) ...`

### insert, delete element to `unordered_map`

- [What is the preferred/idiomatic way to insert into a map?](https://stackoverflow.com/questions/4286670/what-is-the-preferred-idiomatic-way-to-insert-into-a-map)
    ```cpp
    std::unordered_map<int, int> mymap;
    mymap[0] = 42;
    mymap.insert(std::map<int, int>::value_type(0, 42));
    mymap.insert(std::pair<int, int>(1, 42));
    mymap.insert(std::make_pair(2, 42));
    mymap.insert({3, 42});  // list initialization syntax
    mymap.emplace(4, 42);
    ```
- [Remove a key from a C++ map](https://stackoverflow.com/questions/10038985/remove-a-key-from-a-c-map)
    ```cpp
    auto it = mymap.find(0);
    mymap.erase(it);
    mymap.erase(mymap.find(1))
    mymap.erase(2);
    ```
- [Insert or update a map](https://stackoverflow.com/questions/21463928/insert-or-update-a-map)
- [C++ insert, cppreference.com](https://en.cppreference.com/w/cpp/container/map/insert)
    ```cpp
    // if ~ else
    if (map.find(key) == map.end()){
        map.insert(std::pair<int, char>(key, value));
    } else {
        map[key] = value;
    }

    // insert
    // if key exists in map -> the value of success will be 0, and insertion doesn't take place
    // if key is not in map -> the value of success will be 1, and insertion takes place
    auto const result = map.insert(std::make_pair(key, value));
    if (not result.second) { result.first->second = value; }

    const auto [it, success] = m.insert({key, value});
	if (!success) { it->second = value; }

    const auto [it, success] = m.insert({key, value});
	if (!success) { it->second += 1 ; }

    // or_insert_assign
    map.insert_or_assign(key, value);

    // try_emplace
    auto [iterator, inserted] = map.try_emplace(key, value);
    if (!inserted) { iterator->second = value; }
    ```

- [What is the difference? vector::at vs. vector::operator[]](https://stackoverflow.com/questions/9376049/vectorat-vs-vectoroperator)

### iterating over `unordered_map`

- [Iterating over unordered_map of vectors](https://stackoverflow.com/questions/62614386/iterating-over-unordered-map-of-vectors)
    ```cpp
    for (auto &it : _dict){
        std::cout << it.first <<std::endl;
    }

    for (const auto& [key, v] : dict) { .. }
    ```

- iterate over `unordered_map`: `for (const auto& [key, value] : map) { ...  }`

### check if key exists in map

- [How to find if a given key exists in a C++ std::map](https://stackoverflow.com/questions/1939953/how-to-find-if-a-given-key-exists-in-a-c-stdmap)
    ```cpp
    if (m.find(key) == m.end()) {...}
    if (.count(key) > )
    ```

---
## Collections: Unordered_set

- modifying Set during iteration (Be very careful!)
    ```cpp
    // BA02A: motif enumeration
    // wrong! unespected result
    for (auto candidate: candidates) {
        if (ith_nbs.find(candidate) == ith_nbs.end()) {
            candidates.erase(*it)
        }
    }

    // fixed! I wasted a day fixing this part.
    // https://stackoverflow.com/questions/20627458/how-to-remove-elements-from-an-stdset-while-iterating-over-it/20627506#20627506
    for (auto it = candidates.cbegin(); it != candidates.cend();) {
        if (ith_nbs.find(*it) == ith_nbs.end()) {
            it = candidates.erase(it);
        }
        else {
            it++;
        }
    }
    ```

    ```cpp
    #include <iostream>
    #include <vector>
    #include <set>
    using namespace std;

    void vprint(vector<string> vec) {
        for (auto e: vec) {
            cout << e << " ";
        }
        cout << endl;
    }

    void sprint(set<string> myset) {
        for (auto e: myset) {
            cout << e << " ";
        }
        cout << endl;
    }

    int main() {
        // works fine
        vector<string> vec = {"a", "b", "c", "d"};
        for (auto e: vec) {
            vprint(vec);
            vec.pop_back();
        }
        // error: segment failure, sometimes unexpected result
        // hard to find which part is wrong
        set<string> myset = {"a", "b", "c", "d"};
        for (auto e: myset) {  // set iteration
            sprint(myset);
            myset.erase(e);    // modifying set
        }
        // fixed
        // iterator erase( const_iterator first, const_iterator last );  --> 'erase' return iterator
        // https://en.cppreference.com/w/cpp/container/vector/erase
        // TODO: How can I implement a function returning iterator?
        set<string> myset = {"a", "b", "c", "d"};
        for (auto it = myset.begin(); it != myset.end();) {  // set iteration
            sprint(myset);
            it = myset.erase(it++);
        }
        return 0;
    }
    ```

---
## IO

### Read text file

- read text file, store data to vector
    ```cpp
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
    ```
- redirect stdout to file
  - [How do I output all my cout’s to a text file in C++?](https://www.quora.com/How-do-I-output-all-my-cout-s-to-a-text-file-in-C)
    ```cpp
    // print result to text file
    ofstream ofs{ "/home/wsl/rosalind/data/ba03d_output.txt" }; // std::basic_ofstream, <fstream>
    auto backup = cout.rdbuf();                                 // save backup
    cout.rdbuf(ofs.rdbuf());                                    // redirect stdout to ofstream
    { ... }                                                     // do something, cout << ...
    cout.rdbuf(backup);                                         // recover stdout to backup
    ```

### file access

- [How to read input with c++ using program.exe < text.txt](https://stackoverflow.com/questions/10359695/how-to-read-input-with-c-using-program-exe-text-txt)
- [Opening .txt file via executable file (compiled c++ code) on Mac](https://stackoverflow.com/questions/53130998/opening-txt-file-via-executable-file-compiled-c-code-on-mac)
- [embedding a text file in an exe which can be accessed using fopen](https://stackoverflow.com/questions/7366391/embedding-a-text-file-in-an-exe-which-can-be-accessed-using-fopen)
- [run cmake with a txt file](https://stackoverflow.com/questions/69877468/run-cmake-with-a-txt-file)
- [How to set cmake, in order to add txt files into working directory as resource?](https://stackoverflow.com/questions/46995733/how-to-set-cmake-in-order-to-add-txt-files-into-working-directory-as-resource)

---
## Project Structure

### Simple mono-file project without project

- relatively simple, easy
- but... unorganized and duplicated
- how to compile, run:
  - Vscode + Plugin (Code Runner): `"cpp": "cd $dir && g++ $fileName -o $fileNameWithoutExt && $dir$fileNameWithoutExt && rm $dir$fileNameWithoutExt",`
  -  shell command: `g++ -o <NAME> source.cpp`, then `./<NAME>`

### seperation `*.h` and `*.cpp` file, why?

  - [Is it a good practice to place C++ definitions in header files?](https://stackoverflow.com/questions/583255/is-it-a-good-practice-to-place-c-definitions-in-header-files)
  - [Separating class code into a header and cpp file](https://stackoverflow.com/questions/9579930/separating-class-code-into-a-header-and-cpp-file)
  - [Why have header files and .cpp files? [closed]](https://stackoverflow.com/questions/333889/why-have-header-files-and-cpp-files)
  - [What should go into an .h file?](https://stackoverflow.com/questions/1945846/what-should-go-into-an-h-file)

### Multi-file project

#### > Project by Visual Studio

- simple, organized, well structured
- but... IDE is too heavy (>10G)

#### > Project by CMake

- tutorials
  - [C++ Development with Visual Studio Code (C++ extension, CMake Tools extension)](https://www.youtube.com/watch?v=_dXKSi70vJ0)
- structure
  - need `CMakeLists.txt` file
  - [CMake Made Easy: File structure](https://www.ravbug.com/tutorials/cmake-easy/)
  - [C++ Project Structure and Cross-Platform Build With CMake](https://medium.com/swlh/c-project-structure-for-cmake-67d60135f6f5)
- clean (delete all) bulid result
  - [Looking for a 'cmake clean' command to clear up CMake output](https://stackoverflow.com/questions/9680420/looking-for-a-cmake-clean-command-to-clear-up-cmake-output)
- subprojects

---
## External Libraries

### Boost library

- [How to Install Boost Library in C++ on Linux?](https://www.geeksforgeeks.org/how-to-install-boost-library-in-cpp-on-linux/):  `sudo apt-get install libboost-all-dev`

### vcpkg

- [Get started with vcpkg](https://vcpkg.io/en/getting-started.html)
- [Simple C++ project with CMAKE and VCPKG](https://www.youtube.com/watch?v=4z2jmDr36Fc)
