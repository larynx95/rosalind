# Visual Studio Code settings (C++, Rust)

## ▣ Install Languages, Code Runner Plugin
- install vscode in Windows 11
- install plugins: Remote WSL, Code Runner, Go, C/C++, ...
- Code Runner Executor map
  ```
  "cpp": "cd $dir && g++ $fileName -std=c++2a -o $fileNameWithoutExt && $dir$fileNameWithoutExt && rm $dir$fileNameWithoutExt",
  "rust": "cd $dir && cargo fmt && clear && cargo run && cargo clean",
  "javascript": "cd $dir && node $fileName",
  ```

## ▣ Settings for C/C++
- debug or compile, show result, then delete
- `.vscode/c_cpp_properies.json`
  ```
  {
      "configurations": [
          {
              "name": "Linux",
              "includePath": [
                  "${workspaceFolder}/**"
              ],
              "defines": [],
              "compilerPath": "/usr/bin/gcc",
              "intelliSenseMode": "linux-gcc-x64"
          }
      ],
      "version": 4
  }
  ```
- `.vscode/launch.json`
  ```
  {
      // Use IntelliSense to learn about possible attributes.
      // Hover to view descriptions of existing attributes.
      // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
      "version": "0.2.0",
      "configurations": [
          {
              "name": "g++ - Build and debug active file",
              "type": "cppdbg",
              "request": "launch",
              "program": "${fileDirname}/${fileBasenameNoExtension}",
              "args": [],
              "stopAtEntry": false,
              "cwd": "${fileDirname}",
              "environment": [],
              "externalConsole": false,
              "MIMode": "gdb",
              "setupCommands": [
                  {
                      "description": "Enable pretty-printing for gdb",
                      "text": "-enable-pretty-printing",
                      "ignoreFailures": true
                  },
                  {
                      "description": "Set Disassembly Flavor to Intel",
                      "text": "-gdb-set disassembly-flavor intel",
                      "ignoreFailures": true
                  }
              ],
              "preLaunchTask": "C/C++: g++ build active file",
              "miDebuggerPath": "/usr/bin/gdb",
              "postDebugTask": "C/C++: g++ clean",
          }
      ]
  }
  ```
- `.vscode/settings.json`
  ```
  {
      "editor.tabSize": 4,
      "editor.insertSpaces": true,
      "editor.detectIndentation": false,
      "files.associations": { ... },
      "C_Cpp.clang_format_fallbackStyle": "google",
      "C_Cpp.default.cppStandard": "gnu++20",
      "C_Cpp.default.cStandard": "gnu17",
  }
  ```
- `.vscode/tasks.json`
  ```
  {
      "tasks": [
          {
              "type": "cppbuild",
              "label": "C/C++: g++ build active file",
              "command": "/usr/bin/g++",
              "args": [
                  "-fdiagnostics-color=always",
                  "-g",
                  "${file}",
                  "-o",
                  "${fileDirname}/${fileBasenameNoExtension}"
              ],
              "options": {
                  "cwd": "${fileDirname}"
              },
              "problemMatcher": [
                  "$gcc"
              ],
              "group": {
                  "kind": "build",
                  "isDefault": true
              },
              "detail": "Task generated by Debugger."
          },{
              "type": "shell",
              "label": "C/C++: g++ clean",
              "command": "rm",
              "args": [
                "${fileDirname}/${fileBasenameNoExtension}",
              ],
              "presentation": {
                "reveal": "silent",
                "clear": true
              }
            },
      ],
      "version": "2.0.0"
  }
  ```