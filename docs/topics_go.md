# Topics - Go
- Go is simple, fast, easy, convenient, and well supported!

- [Topics - Go](#topics---go)
  - [Numeric data](#numeric-data)
  - [`string` type](#string-type)
    - [No character type in Go](#no-character-type-in-go)
    - [Iteration over characters in a string](#iteration-over-characters-in-a-string)
    - [Substring](#substring)
  - [Slices](#slices)
  - [Iteration, Loop](#iteration-loop)
  - [Map](#map)
  - [Functions](#functions)
    - [Function Argument: pass-by-value vs. pass-by-reference](#function-argument-pass-by-value-vs-pass-by-reference)
    - [Anonymous Function](#anonymous-function)
  - [IO](#io)

## Numeric data

- infinity: `math.Inf(-1)`, `math.Inf(1)`
- random integer
    ```go
    rand.Seed(time.Now().UnixNano())
    mft.Println(rand.Intn(10))  // [0, 10)
    ```

---
## `string` type

### No character type in Go

- Go strings aren't null terminated.
- no character type in Go: byte (uint8), rune (int32), string
- conversion functions: `byte()`, `[]byte`, `rune()`, `string()`

### Iteration over characters in a string

- [Why is rune in golang an alias for int32 and not uint32?](https://stackoverflow.com/questions/24714665/why-is-rune-in-golang-an-alias-for-int32-and-not-uint32)
- The type of `i` in for-loop is rune, not character.

    ```go
    package main

    import (
    	"fmt"
    	"reflect"
    )

    func main() {
    	s := "Hello, 世界"
    	for _, e := range s {
    		fmt.Printf("typeof: %v\n", reflect.TypeOf(e))  // int32 == rune
            fmt.Println(e == rune('H'))                    // true
    		fmt.Println(byte(e) == byte('H'))              // true, false ...
    		fmt.Println(byte(e), []byte(string(e)))        // 76 [231 149 140], ...
    	}
    }
    ```

### Substring

- simple as Python
- [Extracting substrings in Go](https://stackoverflow.com/questions/12311033/extracting-substrings-in-go)
- index-based slice: `str[ <start> : <end> ]` (end, not inclusive)

---
## Slices

- [concatenate two or more slices in Go (merging only)](https://freshman.tech/snippets/go/concatenate-slices/)
    ```go
    s1 := []string{"a", "b", "c"}
    s2 := []string{"x", "y", "z"}
    s3 := []string{"q", "r", "s"}
    s1 = append(s1, s2...) // s1 will be modified
    s1 = append(s1, s3...) // s1 will be modified, again
    ```
- concatenate two or more slices in Go (modifying and merging together): Be very careful!
    ```go
    // wrong way
    rest := cycle[overlap_idx+1: ]   // shallow copied, 'rest' slice points to original 'cycle'
    cycle = append(cycle[:overlap_idx], new_cycle...)
    cycle = append(cycle, rest...)   // <-- 'rest' in this line is not the same as the first 'rest'
    // 'rest' slice points to original 'cycle'

    // correct way
    var rest []string
    rest = append(rest, cycle[overlap_idx+1:]...)
    cycle = append(cycle[:overlap_idx], new_cycle...)
    cycle = append(cycle, rest...)
    ```
- length, capacity of slice (very important)
- [concatenate two slices in Go](https://stackoverflow.com/questions/16248241/concatenate-two-slices-in-go)
    ```go
    arr := []string{"a", "b", "c", "d", "e", "f"}
    res := arr[3:]
    res = append(res, arr[:3]...)
    ```

---
## Iteration, Loop

- classic for-loop: `for i := 0; i < len(text); i++ {...}`
- no `while` in Go: use `for`
    ```go
    // same as while
    i := 0
    for i < 5 {
        fmt.Println("i =", i)
        i++
    }

    // infinite loop
    for { ... }
    ```
- using index, value pair: `for i, e := range <iterator>`
- Variable for value can be omitted. `for i := range <iterator>`

---
## Map

- insert an entry to map
    ```go
    mymap := make(map[string][]string)
    mymap["a"] = []string{"x","y","z"}  // '=' & ':='
    ```
- [iterating over all the keys of a map](https://stackoverflow.com/questions/1841443/iterating-over-all-the-keys-of-a-map)
    ```go
    for k, v := range m {
        fmt.Printf("key[%s] value[%s]\n", k, v)
    }
    ```

- get all keys, values from a map
  - [Getting a slice of keys from a map](https://stackoverflow.com/questions/21362950/getting-a-slice-of-keys-from-a-map)

    ```go
    keys := make([]int, len(mymap))
    i := 0
    for k := range mymap {
        keys[i] = k
        i++
    }
    ```

- check if key in map
  - [How to check if a map contains a key in Go?](https://stackoverflow.com/questions/2050391/how-to-check-if-a-map-contains-a-key-in-go)
    ```go

    v, ok := map[key]         // map[key] returns (value, true) if key is in map, else returns (default value, false)
    if !ok { ... }            // if absent

    if valuei, ok := pairMap[v]; ok { ... }   // if key is in map
    if valuei, ok := pairMap[v]; !ok { ... }  // if key is not in map

    // map[key] is not value in Go (BA03D, BA03E)
    if values, ok := graph[prefix]; ok {
        values = append(values, suffix)
        // graph[prefix] = append(graph[prefix], suffix)  // error
    } else {
        graph[prefix] = []string{suffix}
    }
    ```

---
## Functions

### Function Argument: pass-by-value vs. pass-by-reference

- [passing by reference and value in Go to functions](https://stackoverflow.com/questions/47296325/passing-by-reference-and-value-in-go-to-functions)

### Anonymous Function

- [Go by Example: Closures](https://gobyexample.com/closures)


---
## IO

- Redirect `os.Stdout` to file
    ```go
    old := os.Stdout // save os.Stdout
    f, err := os.Create("/home/wsl/rosalind/data/ba03d_output.txt")
    if err != nil {
        log.Fatal(err)
    }
    defer f.Close()
    os.Stdout = old // recover os.Stdout
    ```
