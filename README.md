This is the original version of Sapphire by [Mary Beard](https://isnap.nd.edu/people/group-pages/mary-beard/). 
Only the README-file has been changed.


- [General Information](#general-information)
- [Compilation & Execution](#compilation--execution)
  - [Requirements](#requirements)
  - [Build](#build)
  - [Run](#run)

--- 
# General Information

# Compilation & Execution
## Requirements
Compiling Sapphire requires `CMake`, `ROOT`, and `GNU Scientific Library (GSL)`. Check the files in `cmake/Modules` whether cmake is looking for these packages in the right place. 

For an `automized generation of documentation files`, one needs `doxygen` installed on the system.

For the `Message Passing Interface (MPI)` build, `boost 1.4` libraries are required. The MPI build is optional and can be enabled in the `CMakeLists.txt`.

## Build

1. Create a build directory under the main Sapphire directory and change to that directory. 
   ```bash
   mkdir build
   cd build
   ```
2. Run CMake against the main directory, optionally specifying the desired C++ compiler.
    ```bash
    cmake -DCMAKE_CXX_COMPILER=icpc  -DCMAKE_C_COMPILER=icc 
    ```
3. Type make install to build Sapphire.  The executable is put in
  the build directory.
    ```bash
    make install
    ```

Be aware that Sapphire links the paths to the needed tables at compile time.
While the executable can be moved, the main Sapphire directory should stay in
place.  If moved, the build process should be repeated.

## Run

To execute the code just enter

```
./sapphire X+a
```

where X is the heavy nucleus and a is the projectile. Examples are 25Mg+a or 60Fe+n.  
