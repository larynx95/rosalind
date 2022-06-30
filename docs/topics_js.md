# Topics - JavaScript

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
