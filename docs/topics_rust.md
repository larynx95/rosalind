# Topics - Rust
- Don't waste time to learn Rust!
- If Rust can, other languages can do it. Even easier and faster.

- [Topics - Rust](#topics---rust)
  - [Numbers](#numbers)
  - [Strings](#strings)
    - [`String`, `&String`, `str`, `&str`](#string-string-str-str)
    - [`char`, `byte`](#char-byte)
    - [string join](#string-join)
    - [string split](#string-split)
    - [string length, count](#string-length-count)
    - [character to string](#character-to-string)
  - [Iteration](#iteration)
    - [Modifying iterator during looping](#modifying-iterator-during-looping)
    - [What is the type of local variable `i` in `for-loop`](#what-is-the-type-of-local-variable-i-in-for-loop)
    - [`iter()`, `into_iter()`, `inter_mut()`](#iter-into_iter-inter_mut)
    - [while-let](#while-let)
  - [Data Structure: Array](#data-structure-array)
  - [Data Structure: Tuple](#data-structure-tuple)
  - [Data Structure: Vector](#data-structure-vector)
    - [is_empty](#is_empty)
    - [remove element](#remove-element)
    - [get subvector](#get-subvector)
    - [mutaing vector](#mutaing-vector)
    - [merge two or more vectors](#merge-two-or-more-vectors)
    - [iterating over vector](#iterating-over-vector)
  - [Data Structure: HashMap](#data-structure-hashmap)
    - [insert/remove element in HashMap](#insertremove-element-in-hashmap)
    - [check if key exists in HashMap](#check-if-key-exists-in-hashmap)
    - [updating values of HashMap](#updating-values-of-hashmap)
    - [iterating over HashMap](#iterating-over-hashmap)
  - [Data Structure: HashSet](#data-structure-hashset)
  - [IO](#io)
  - [Ownership, Borrowing, Lifetime](#ownership-borrowing-lifetime)
    - [Mutable](#mutable)
    - [Ownership](#ownership)
    - [Borrowing](#borrowing)
    - [Lifetime](#lifetime)
    - [Function Argument: call-by-value vs. call-by-reference](#function-argument-call-by-value-vs-call-by-reference)
  - [Wrapper: Cell, RefCell, Rc, Arc, Mutex, Box, RwLock, Cow](#wrapper-cell-refcell-rc-arc-mutex-box-rwlock-cow)
  - [Structure, Enum, Impl, Trait](#structure-enum-impl-trait)
  - [Project, Module, pub, use ...](#project-module-pub-use-)
  - [Recursion](#recursion)

## Numbers

- pow function: `u64pow(u32)`
- only `usize` in for-loop
- size of vector, array: `let freq_vec = vec![0; 4u32.pow(k as u32) as usize];`
- infinity: `std::f64::INFINITY`, `std::usize::INFINITY`
- random integer
    ```rs
    let mut rng = rand::thread_rng();
    let n1: u8 = rng.gen();
    let n2: u16 = rng.gen();
    // one line
    let num = rand::thread_rng().gen_range(0..100);  // [0..100)
    ```

---
## Strings

### `String`, `&String`, `str`, `&str`

- [What are the differences between Rust's `String` and `str`?](https://stackoverflow.com/questions/24158114/what-are-the-differences-between-rusts-string-and-str)
  - `String` is the dynamic heap string type, like Vec: use it when you need to own or modify your string data.
  - `str` is an immutable1 sequence of UTF-8 bytes of dynamic length somewhere in memory. Since the size is unknown, one can only handle it behind a pointer. This means that str most commonly2 appears as &str: a reference to some UTF-8 data, normally called a "string slice" or just a "slice". A slice is just a view onto some data, and that data can be anywhere, e.g.

### `char`, `byte`

- `char` type exists in Rust. (UTF)
- many kinds of string types: string, &string, str, &str
- many kinds of character types: byte, char
- [How to get the last character of a &str?](https://stackoverflow.com/questions/48642342/how-to-get-the-last-character-of-a-str)
    ```rs
    mystring.chars().last().unwrap();
    ```

### string [join](https://docs.rs/itertools/0.6.0/itertools/fn.join.html)

- [join](https://docs.rs/itertools/0.6.0/itertools/fn.join.html)
- [What is the equivalent of the join operator over a vector of Strings?](https://stackoverflow.com/questions/28311868/what-is-the-equivalent-of-the-join-operator-over-a-vector-of-strings)
    ```rs
    let string_list = vec!["Foo".to_string(),"Bar".to_string()];
    let joined = string_list.join("-");
    assert_eq!("Foo-Bar", joined);
    ```

### string split

- [How do I split a string in Rust?](https://stackoverflow.com/questions/26643688/how-do-i-split-a-string-in-rust)
    ```rs
    let mut split = "some string 123 ffd".split("123");
    for s in split {
        println!("{}", s)
    }
    let vec = split.collect::<Vec<&str>>();
    // OR
    let vec: Vec<&str> = split.collect();
    ```

### string length, [count](https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.count)

- all methods of [primitive::str](https://doc.rust-lang.org/std/primitive.str.html)
- [Get the String length in characters in Rust](https://stackoverflow.com/questions/46290655/get-the-string-length-in-characters-in-rust)
    ```rs
    println!("{}", "string".chars().count());
    ```

### character to string

- [How to concatenate a char onto a string in Rust?](https://stackoverflow.com/questions/37889337/how-to-concatenate-a-char-onto-a-string-in-rust)
    ```rs
    // push
    let mut a_string = String::from("Hello World");
    a_string.push('!');

    // format
    let s = String::from("March");
    let s1 = format!("{}!", s);         // example 1
    let s2 = format!("{}{}", s, '!');   // example 2
    println!("{} {}", s1, s2);          // print
    ```

---
## Iteration

### Modifying iterator during looping

### What is the type of local variable `i` in `for-loop`

- [What is the purpose of `&` before the loop variable?](https://stackoverflow.com/questions/57339201/what-is-the-purpose-of-before-the-loop-variable)
- The local variable `i` is `&i32`.

### `iter()`, `into_iter()`, `inter_mut()`

- `iter()`: borrow, take reference of each elem, original collection untouched
- `into_iter()`: consume, move ownership, can't use it again
- `iter_mut()`: mutably borrow

### while-let

---
## Data Structure: Array

- Is Array really necessary in Rust?
- [Create empty array with a length from a variable [duplicate]](https://stackoverflow.com/questions/44847574/create-empty-array-with-a-length-from-a-variable)
- [Is it possible to have stack allocated arrays with the size determined at runtime in Rust?](https://stackoverflow.com/questions/27859822/is-it-possible-to-have-stack-allocated-arrays-with-the-size-determined-at-runtim)

---
## Data Structure: Tuple

- insert/remove/change tuple
    ```rs
    // HashMap
    let mut mymap: HashMap<String, (usize, usize)> = HashMap::from([
        ("a".to_string(), (1,1)),
        ("b".to_string(), (2,2)),
        ("c".to_string(), (3,3)),
    ]);
    // declare, initialize Tuple
    let t1: (usize, usize) = (4,4);
    // insert
    mymap.insert("d".to_string(), t1);
    mymap.insert("e".to_string(), (5,5));
    // chamge value of tuple
    let mut t2: (usize, usize) = (5,5);
    t2.1 += 100;
    mymap.insert("f".to_string(), t2);
    println!("{:?}", mymap);
    ```

---
## Data Structure: Vector

### is_empty

- [is_empty](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.is_empty)
    ```rs
    if vec.is_empty() {...}
    if vec.len() == 0 {...}
    ```

### remove element

- insert, remove element from Vec by index, value
- remove element at an index: [Remove an element from a vector](https://stackoverflow.com/questions/26243025/remove-an-element-from-a-vector)
    ```rs
    // insert element to Vec
    vec.push(elem)
    // insert element to Vec, at an index
    vec.insert(index, elem)
    // remove element from Vec
    let removed_elem = vec.remove(index)
    // inter -> position
    let index = vec.iter().position(|x| *x == some_x).unwrap();
    vec.remove(index);
    ```
- remove element by value
- remove & retrieve: [pub fn pop(&mut self) -> Option<T>](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.pop)

### get subvector

- [How do I get a slice of a Vec<T> in Rust?](https://stackoverflow.com/questions/39785597/how-do-i-get-a-slice-of-a-vect-in-rust)
    ```rs
    let a = vec![1, 2, 3, 4, 5];
    // With a start and an end
    println!("{:?}", &a[1..4]);  // &[{integer}]
    println!("{:?}", a[1..4].to_vec());  // Vec<int>
    ```
- filtering
  ```rs
  let v = vec!["aaa", "bbb", "ccc", "ddd"];
  let f = v.iter().enumerate().filter(|&(i, e)| i != 0).map(|(_,e)| *e).collect::<Vec<&str>>();
  println!("{:?}", f);
  ```

### mutaing vector

- [retain](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.insert), filter
    ```rs
    // retain
    let mut vec = vec![1, 2, 3, 4];
    vec.retain(|&x| x % 2 == 0);
    println!("{:?}", vec);  // [2,4]

    // filter
    let mut vec = vec![1,2,3,4];
    vec.iter().filter(|&x| x%2 == 0);
    println!("{:?}", vec);  // [1,2,3,4]

    // vec -> iter -> filter -> collect
    let mut vec = vec![1,2,3,4];
    vec = vec.iter().filter(|&x| x%2 == 0).map(|x| *x).collect::<Vec<usize>>();
    println!("{:?}", vec);  // [2,4]
    ```

### merge two or more vectors

- merge two vectors ([Best way to concatenate vectors in Rust](https://stackoverflow.com/questions/40792801/best-way-to-concatenate-vectors-in-rust))
    ```rs
    // mutating both vectors
    let mut a = vec![1, 2, 3];
    let mut b = vec![4, 5, 6];
    a.append(&mut b);
    assert_eq!(a, [1, 2, 3, 4, 5, 6]);
    assert_eq!(b, []);

    // mutating a, move b
    let mut a = vec![1, 2, 3];
    let b = vec![4, 5, 6];
    a.extend(b);
    assert_eq!(a, [1, 2, 3, 4, 5, 6]);
    // b is moved and can't be used anymore

    // mutating a, b intact
    let mut a = vec![1, 2, 3];
    let b = vec![4, 5, 6];
    a.extend(&b);
    assert_eq!(a, [1, 2, 3, 4, 5, 6]);
    assert_eq!(b, [4, 5, 6]);
    ```

### iterating over vector

- [How to iterate a Vec<T> with the indexed position? (like Golang)](https://stackoverflow.com/questions/28991050/how-to-iterate-a-vect-with-the-indexed-position)
    ```rs
    fn main() {
        let v = vec![1; 10];
        for (pos, e) in v.iter().enumerate() {
            println!("Element at position {}: {:?}", pos, e);
        }
    }
    ```

---
## Data Structure: HashMap

### insert/remove element in HashMap

- remove: [`retain`](https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.retain), [`remove`](https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.remove), [`remove_entry`](https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.remove_entry)
    ```rs
    // insert
    map.insert(key, value);
    // insert, remove
    let mut map = HashMap::from([(1,"a"), (3,"c"), (4,"d")]);
    map.insert(2,"b");                        // insert
    map.retain(|&k, _| k != 2);               // remove by retain
    let removed_key = map.remove(&1);         // remove & return key
    let removed_entry = map.remove_entry(&3); // remove & return entry
    println!("{:?}", removed_key);            // Some("a")
    println!("{:?}", removed_entry);          // Some((3,"c"))
    ```

### check if key exists in HashMap

- [How to lookup from and insert into a HashMap efficiently?](https://stackoverflow.com/questions/28512394/how-to-lookup-from-and-insert-into-a-hashmap-efficiently)
- [pub fn entry(&mut self, key: K) -> Entry<'_, K, V>](https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.entry)
    ```rs
    map.contains_key(key)  // check if key is in hashMap
    // entry
    let values = match map.entry(key) {
        Entry::Occupied(o) => o.into_mut(),
        Entry::Vacant(v) => v.insert(default),
    };
    let values = map.entry(key).or_insert_with(|| default);  // or_insert_wwith
    let values = map.entry(key).or_insert(default);          // or_insert
    let values = map.entry(key).or_default();                // or_default
    ```

### updating values of HashMap

- [StackOverflow: How can I update a value in a mutable HashMap?](https://stackoverflow.com/questions/30414424/how-can-i-update-a-value-in-a-mutable-hashmap)
    ```rs
    *my_map.get_mut("a").unwrap() += 10;
    *my_map.entry("a").or_insert(42) += 10;
    my_map.inert(key, my_map[key] + 10);
    ```
- modifying vector value of HashMap (BA02D)
    ```rs
    // change value (vector)
    let mut map = HashMap::from([('a', vec![1,2,3]),('b', vec![4,5,6]),('c', vec![7,8,9])]);
    println!("{:?}", map.get_mut(&'a'));
    println!("{:?}", map.get_mut(&'t'));
    map.get_mut(&'a').unwrap().push(100);  // get_mut()
    println!("{:?}", map);
    // change type of value - 220613
    let mut map = HashMap::from([('a', vec![1,2,3]),('b', vec![4,5,6]),('c', vec![7,8,9])]);
    for (key, val) in map.iter() {
        // get key, value
        println!("{:?}", map.get(&key).unwrap());
        println!("{:?}", val);
        //
        let vec_float: Vec<f64> = val.iter().map(|&x| x as f64 / 10.0).collect::<Vec<_>>();
        println!("{:?}", vec_float);
    }
    ```

### iterating over HashMap

- (key-value) pair, key
    ```rs
    // iter()
    for (key, value) in map.iter() {...}  // Hashmap can be used again later.
    // just HashMap
    for (key, value) in map {...}         // HashMap moved here!
    ```

---
## Data Structure: HashSet

- Intersection of HashSet
  - [How do I intersect two HashSets while moving values in common into a new set?](https://stackoverflow.com/questions/55975234/how-do-i-intersect-two-hashsets-while-moving-values-in-common-into-a-new-set)
  - [Struct std::collections::HashSet::intersection](https://doc.rust-lang.org/std/collections/struct.HashSet.html#method.intersection)
  - [retain: pub fn retain<F>(&mut self, f: F)](https://doc.rust-lang.org/std/collections/struct.HashSet.html#method.retain)
  - [pub fn intersection<'a>](https://doc.rust-lang.org/std/collections/struct.HashSet.html#method.intersection)
    ```rs
    use std::collections::HashSet;
    let a = HashSet::from([1, 2, 3]);
    let b = HashSet::from([4, 2, 3, 4]);
    // Print 2, 3 in arbitrary order.
    for x in a.intersection(&b) {
        println!("{x}");
    }
    let intersection: HashSet<_> = a.intersection(&b).collect();
    assert_eq!(intersection, [2, 3].iter().collect());
    ```
- [How to join elements of HashSet into a String with a delimiter](https://stackoverflow.com/questions/47578011/how-to-join-elements-of-hashset-into-a-string-with-a-delimiter)
    ```rs
    println!("{}", hash_set.iter().join(", "));
    println!("{}", itertools::join(&hash_set, ", "));

    // BA03D
    for (key, values) in debruijn_from_text(text, k) {
        print!("{key} -> ");
        print!("{}", values.into_iter().collect::<Vec<String>>().join(","));
        println!()
    }
    ```
- delete a element in a value
    ```rs
    let mut map = HashMap::from([
        ("a", 1),
        ("b", 2),
        ("c", 3),
    ]);
    map = map.iter().filter(|&(k,v)| k != &"a").map(|(k,v)| (*k, *v)).collect::<HashMap<&str, usize>>();
    println!("{:?}", map);  // {"c": 3, "b": 2}
    ```

---
## IO

- read text file
  ```rs
  // helper function: read a text file, create a Vector of strings
  fn lines_from_file(filename: &str) -> Vec<String> {
      let file = File::open(filename).expect("no such file");
      let buf = BufReader::new(file);
      buf.lines()
          .map(|l| l.expect("Could not parse line"))
          .collect()
  }
  ```
- write to text file
  ```rs
  // write to text file, example (BA01N)
  let mut w = File::create("output.txt").unwrap();
  let neighborhood = iterative_neighbors(pattern, k);
  for nb in neighborhood {
      writeln!(&mut w, "{}", nb).unwrap();
  }
  ```

---
## Ownership, Borrowing, Lifetime

### Mutable

- [What's the difference between placing "mut" before a variable name and after the ":"?](https://stackoverflow.com/questions/28587698/whats-the-difference-between-placing-mut-before-a-variable-name-and-after-the)
    ```
    // Rust          C/C++
        a: &T     == const T* const a; // can't mutate either
    mut a: &T     == const T* a;       // can't mutate what is pointed to
        a: &mut T == T* const a;       // can't mutate pointer
    mut a: &mut T == T* a;             // can mutate both
    ```

### Ownership
- Resource acquisition is initialization (RAII) in C++
- [The Rust Programming Language: Understanding Ownership](https://doc.rust-lang.org/book/ch04-00-understanding-ownership.html)

### Borrowing

### Lifetime

### Function Argument: call-by-value vs. call-by-reference

- [Why does Rust have both call by value and call by reference?](https://stackoverflow.com/questions/36562262/why-does-rust-have-both-call-by-value-and-call-by-reference)
- [Is using `ref` in a function argument the same as automatically taking a reference?](https://stackoverflow.com/questions/21575198/is-using-ref-in-a-function-argument-the-same-as-automatically-taking-a-referen)

---
## Wrapper: Cell, RefCell, Rc, Arc, Mutex, Box, RwLock, Cow
## Structure, Enum, Impl, Trait
## Project, Module, pub, use ...
## Recursion

- [Why is recursion not suggested in Rust?](https://stackoverflow.com/questions/65948553/why-is-recursion-not-suggested-in-rust)
