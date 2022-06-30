# Comparison
my summary

## String

### Get substring from a string
  | lang       | expression              |
  | ---------- | ----------------------- |
  | Python     | `str[i : i+k]`          |
  | Go         | `str[i : i+k]`          |
  | JavaScript | `str.substring(i, i+k)` |
  | C++        | `str.substr(i, k)`      |
  | Rust       | `&str[i..(i+k)]`        |

### Chracter type
  | lang       | expression                                                 |
  | ---------- | ---------------------------------------------------------- |
  | Python     | no char type, `<class 'str'>`, only string with length one |
  | Go         | no char type, `byte`(ASCII) or `rune`(UTF)                 |
  | JavaScript | no char type, `string`                                     |
  | C++        | `Char`(one byte == 8 bits)                                 |
  | Rust       | `char`(UTF)                                                |

### Chracter == Integer

---
## Numbers

### Infinity
  | lang       | expression                                              |
  | ---------- | ------------------------------------------------------- |
  | Python     | `float('inf')`, `float(-inf')`, `math.inf`, `-math.inf` |
  | Go         | `math.Inf(-1)`, `math.Inf(1)`                           |
  | JavaScript | `Infinity`                                              |
  | C++        | `std::numeric_limits<double>::infinity()`               |
  | Rust       | `usize::MAX`, `usize::MIN`                              |

- [Setting an int to Infinity in C++](https://stackoverflow.com/questions/8690567/setting-an-int-to-infinity-in-c)
  - `std::numeric_limits<int>::max()`
  - `std::numeric_limits<double>::infinity()`
- [Rust-lang.org: usize](https://doc.rust-lang.org/std/primitive.usize.html)

### Default numeric types

---
## Loop, Iterator, Generator

### Loops
  | lang       | expression                                                  |
  | ---------- | ----------------------------------------------------------- |
  | Python     | `for`, `for-range`, `while` + comprehension + HOF           |
  | Go         | `for`, `for-range`, `while` (simple is best)                |
  | JavaScript | `for`,`while` + `forEach` + HOF + cahining                  |
  | C++        | `for`, `for-range`, while + HOF(transform,accumulate        |
  | Rust       | `for-range`, `loop`, `while` + iterators + `foreach` +  HOF |

---
## Functions

### Arument passing: value, reference
  | lang       | expression               |
  | ---------- | ------------------------ |
  | Python     | by object reference      |
  | Go         | by value                 |
  | Java       | by value                 |
  | JavaScript | by value                 |
  | C          | by value                 |
  | C++        | both value and reference |
  | Rust       | by value                 |

- ["All arguments to functions are passed by value" in C, confusion about pass by reference in C++](https://stackoverflow.com/questions/41413124/all-arguments-to-functions-are-passed-by-value-in-c-confusion-about-pass-by-r)
- [Understanding JavaScript Pass-By-Value](https://www.javascripttutorial.net/javascript-pass-by-value/)
- [How are arguments passed by value or by reference in Python?](https://www.tutorialspoint.com/how-are-arguments-passed-by-value-or-by-reference-in-python)

---
## Vector Operations

### Get subvectors

  | lang       | expression                                          |
  | ---------- | --------------------------------------------------- |
  | Python     |                                                     |
  | Go         |                                                     |
  | JavaScript | `ar.slice(1, ar.length-1)`                          |
  | C++        | `vector<T> subvec { vec.begin() + idx, vec.end() }` |
  | Rust       | `vec[1..vec.len()-1].to_vec()`                      |

---
## Map Operations

### Check if key exists in a map
  | lang       | expression                                              |
  | ---------- | ------------------------------------------------------- |
  | Python     | `if key in dictionary:`                                 |
  | Go         | `if val, ok := dict["foo"]; ok { ... }`                 |
  | JavaScript | `if (map.has(key))`                                     |
  | C++        | `if (map.count(key) > 0)`, `if (m.find(ch) != m.end())` |
  | Rust       | `if map.contains_key(&key)`                             |

- [Rust: How to lookup from and insert into a HashMap efficiently?](https://stackoverflow.com/questions/28512394/how-to-lookup-from-and-insert-into-a-hashmap-efficiently)

### Get value by key

### Iterate over Map
  | lang       | expression                                                   |
  | ---------- | ------------------------------------------------------------ |
  | Python     |                                                              |
  | Go         |                                                              |
  | JavaScript |                                                              |
  | C++        | `for (const auto& [key, v] : map)`, `for (auto &it : _dict)` |
  | Rust       |                                                              |

---
## Set Operations
  | lang       | expression |
  | ---------- | ---------- |
  | Python     |            |
  | Go         |            |
  | JavaScript |            |
  | C++        |            |
  | Rust       |            |