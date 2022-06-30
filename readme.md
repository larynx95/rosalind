# Rosalind
Rosalind (2022.5) - Python, JavaScript, Go, Rust, C++ ... for hobby

- [Rosalind](#rosalind)
  - [Install Programming Languages](#install-programming-languages)
    - [WSL, ZSH](#wsl-zsh)
    - [Languages, Compilers](#languages-compilers)
    - [Best languages](#best-languages)
  - [My Learning Curve](#my-learning-curve)
    - [Common, Comparison](#common-comparison)
    - [Docs, Summary, Language specific](#docs-summary-language-specific)
    - [Progress](#progress)
  - [Usefule Sites](#usefule-sites)
    - [C/C++](#cc)
    - [Go](#go)
    - [Rust](#rust)
    - [Python](#python)
    - [JavaScript](#javascript)

## Install Programming Languages

### WSL, ZSH
- install WSL in Windows 11
- update, upgrade Ubuntu: `sudo apt update -y; sudo apt upgrade -y`
- install ZSH shell, Oh-My-ZSH, and it's plugins

### Languages, Compilers
- **Python**: already installed in WSL Ubuntu
- **nodejs**: [Microsoft: Install Node.js on Windows Subsystem for Linux (WSL2)](https://docs.microsoft.com/en-us/windows/dev-environment/javascript/nodejs-on-wsl)
  - `curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.1/install.sh | bash`
  - `nvm install --lts`
- **Go**: [install Go](https://go.dev/doc/install), download [latest Go version](https://go.dev/dl/go1.18.3.linux-amd64.tar.gz), extract to `/usr/local`
- **Rust**: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
- **C++**: `sudo apt install gcc g++ make cmake gdb -y`
- [VScode settings for C/C++](https://github.com/larynx95/rosalind/blob/main/docs/config_vscode.md)

### Best languages
- Python >>> Javascript > Go >>> C++ >>> Rust

---
## My Learning Curve

### Common, Comparison
- [Common Skills for Rosalind](https://github.com/larynx95/rosalind/blob/main/docs/common.md)
- [Comparison](https://github.com/larynx95/rosalind/blob/main/docs/comparison.md)

### Docs, Summary, Language specific
- [Python](https://github.com/larynx95/rosalind/blob/main/docs/topics_py.md)
- [JavaScript](https://github.com/larynx95/rosalind/blob/main/docs/topics_js.md)
- [Go](https://github.com/larynx95/rosalind/blob/main/docs/topics_go.md)
- [C++](https://github.com/larynx95/rosalind/blob/main/docs/topics_cpp.md)
- [Rust](https://github.com/larynx95/rosalind/blob/main/docs/topics_rust.md)

### Progress

- progress
    | langs      | progress | todo                    |
    | ---------- | -------- | ----------------------- |
    | Python     | ba06     | ba03k, ba04im, ba05jklm |
    | JavaScript | ba04     | ba02g, ba04eim          |
    | Go         | ba03h    | ba02g                   |
    | C++        | ba03h    | ba01n, ba02f            |
    | Rust       | ba03h    | ba01n                   |

- folder structure
    ```
    ┌ .vscode             vscode settings
    ├ data ...            data text file
    ├ docs ...            summary of idioms ... markdown files
    ├ include             C++ header
    │   utils.hpp
    └ src                 source codes
       ├ ba01             chapter
       │   ├ ba01a        each problem
       │   ├ ba01b
       │   ├ ...
       │   └ ba01n
       ├ ba02 ...
       ├ ...
       └ ba11 ...
      .gitignore         git ignore file
      readme.md          this file
      toc_rosalind       table of contents (list of problems)
---
## Usefule Sites

### C/C++
- [Techie Delight, online C/C++ IDE](https://techiedelight.com/compiler/)
- [cdecl: pointer decoder](https://cdecl.org/)
- [C++ references](https://en.cppreference.com/w/)

### Go
- [Go playground](https://go.dev/play/)
- [Go docs](https://pkg.go.dev/google.golang.org/api/docs/v1)

### Rust
- [Rust playground](https://play.rust-lang.org/)
- [Crate std](https://doc.rust-lang.org/stable/std/)
- [Rust by Example](https://doc.rust-lang.org/rust-by-example/)

### Python
- [Python visualizer](https://pythontutor.com/visualize.html)
- [Code skulptor](https://py3.codeskulptor.org/)

### JavaScript
- [Programiz, online JavaScript IDE](https://www.programiz.com/javascript/online-compiler/)