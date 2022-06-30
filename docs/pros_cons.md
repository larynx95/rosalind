# Pros and Cons
What I felt...

## Python for Rosalind

### 👍 Pros
- the BEST for Rosalind, no doubt
- battery included: list/set comprehension, HOF, generator, recursion ...
- multi-purpose, easy to learn

### 👎 Cons
- slow, nothing else

---
## Javascript for Rosalind

### 👍 Pros
- much faster than Python, sometimes faster than compiling languages
- elegant syntax
- suitable for expressing much more diverse code
- a language that evolves very quickly

### 👎 Cons
- hard to find what is wrong in the code
- easy to use, but more difficult to learn than Python
- **many unexpected results**: comparison, equality checking
  - array (object) comparison problem
  - object values or keys are not unique even in Set/Map/Object
  - complex equality checking, even primitive types
  - object property must be string type
  - sorting a integer array by built-in `sort()` function (very weird result!)
- emptiness checking
  - a Set/Map/Object of lists/maps/sets
  - condition for emptiness
- loop
  - `For-in` loops is not what I though of: `for-in`, `for-of`, `for-each` ... very confusing

---
## Go for Rosalind

### 👍 Pros
- faster than Python
- better than C, garbage collector, compiling, strict typed
- garbage collector, but very fast

### 👎 Cons
- pointer, but simpler than C/C++
- confusing string type: string vs. []string
- string handling: character, byte, rune, UTF-8, slice, substring, concatenation
- no Set (but Map can do the same thing)

---
## Rust for Rosalind

### 👍 Pros
- great compiler

### 👎 Cons
- There're much better options. Python, JavaScript, Golang
  - Python: the best of all (for Rosalind)
  - JavaScript: much faster than Python, even faster than Rust, sometimes faster than C++
  - Golang: if speed matters, Go (C/C++)
- Complex, difficult, steep learning curve
- many option for doing the same thing (distraction)
- Without full understanding, coding can be unnecessarily painful and frustrating.
- Even simple things can be complicated and verbose.
- Not so fast as I expected!

---
## C++ for Rosalind

### 👍 Pros
- fast, maybe Fastest, undoubtedly legendary programming language
- many gurus, well supported ecosystem, a language that creates another language
- not old, even modern (recently updated C++17), contains almost everything, used for everything

### 👎 Cons
- poor readability
- difficult, too much to learn, heavy keeping on gaining weight
- complex compiling in multifile, large scale project
- A project with multiple files is very difficult to build without Visual Studio. CMake tool is not easy for beginners.
- modifying iterator duing for-loop problem (including set) - TODO: Compare to other languages.
    ```cpp
    // BA02A: motif enumeration
    // no error! but wrong, unespected result, hard to find what is wrong
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
