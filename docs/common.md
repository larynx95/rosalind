# Language Comparison
Common skills to solve Rosalind quizzes

---
## Types
### Numeric data

- most frequently used numeric types: int, double
- positive/negative infinity, max/min value
- random number

### Constness

- JavaScript: `const`, `let`, `var`
- C++: `const`, `constexpre`

---
## Strings

### String manipulations

- character type: byte, int, rune, ... and comparison between character and integer
- index-based operation (especially in Rust), get character and substring
- length of a string
- iteration over characters in a string
- pair(indices, character), enumerate, zip functions

### Type conversions

- type conversion: string <-> int, float, double
- storing string in stack, heap

### Iteration over string

- classic loops (for, while),
- range-based for, foreach, Higer Order Fucntions(HOF)
- Simple is the best.

---
## Iterations, Looping, Generators, Conditionals

### Classic `for`. `while`

- classic three components for-loop
- range-based for-loop

### Iterator

- iterable, iterator
- (Rust) mutable/move/borrow iterators: `iter()`, `iter_mut()`, `into_iter()`
- (Rust) destructive/moved iterator
- (Rust) `while let`

### Generator

- Python, JavaScript: `yield`
- Go, C++: complex, not easy
- Rust: [Rust unstable book: generators](https://doc.rust-lang.org/beta/unstable-book/language-features/generators.html)

### Recursion, Dynamic Programming

- Python: tail-recursion, recursion depth problem
- (Rust) Recursion is not recommended in Rust. TODO: Why?

### Conditionals

- classic `if-else`, `elseif`
- (Rust) `match`, `if let`, other usefule methods `or_insert`...

---
## Memory Management

### Garbage collector, Ownership/Borrowing, Manual

- C/C++: manual
- JavaScript, Python, Go: Garbage collector
- Rust: Ownership, Borrowing

### Pointers, References

- Python, JavaScript: no pointer
- Go, Rust: has pointer, reference
- C/C++: full of pointers (C/C++), references (C++)

### Call-by-value, Call-by-reference

### Pointers, References

- C: full of pointers
- C++: pointers + references
- Rust: reference + move/borrow/mut...
- Go: pointers, automatic dereferencing

---
## Functions

### Function Arguments, Return value

- mutable/immutable, call-by-value/reference
- JavaScript: All function arguments are always passed by value.
- C: All arguments to functions are passed by value in C.
- C++:

### Anonymous function, Closure

- sometimes very useful

### HOF: Map, Filter, Reduce, Fold

---
## Array

- array size at compile time problem

## Vector, Slice

- length, capacity, original/copied
- create, add, remove, change element
- subvector, part of slice
- shallow/deep copy concepts

---
## HashMap, Dictionary, Map

- declaring HashMap (ordered, unordered), declare and initialize HashMap simultaneously
- reset whole HashMap
- add, remove, change key-value pairs in a HashMap
- comapre two HashMaps
- update/modify value (simple, complex types)
- (JavaScript) Object is not a HashMap. The type of Object's properties is string. Other types are not allowed.

---
## HashSet, Set

- declaring HashSet, declare and initialize HashSet simultaneously
- reset whole HashSet
- add, remove, change key-value pairs in a HashSet
- comapre two HashSets
- HashSet can't be changed during iteration!
- (Go) There's no Set in Golang. Map can do the same thing.

---
## Project Structure

- Go, Rust need separation (go mod, cargo.toml)
- Script languages no not need separation.
- C/C++: Visual studio, CMake
