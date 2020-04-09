This is the original version of Sapphire by [Mary Beard](https://isnap.nd.edu/people/group-pages/mary-beard/). 
Only the README-file has been changed.


- [General Information](#general-information)
- [Compilation & Execution](#compilation--execution)
  - [Requirements](#requirements)
  - [Build](#build)
  - [Run](#run)
- [Features](#features)
- [References](#references)
  - [About Sapphire](#about-sapphire)
  - [Using Sapphire](#using-sapphire)

--- 
# General Information

# Compilation & Execution
## Requirements
Compiling Sapphire requires [CMake](https://cmake.org/), [ROOT](https://root.cern.ch/), and [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). Check the files in `cmake/Modules` whether cmake is looking for these packages in the right place. 

For an `automized generation of documentation files`, one needs [doxygen](http://doxygen.nl/) installed on the system.

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

# Features

## Calculate Cross Sections



# References
## About Sapphire
If you are using *Sapphire* you should cite one or more of the references below.

* M. Beard, E. Uberseder, R. Crowter, and M. Wiescher, *Comparison of statistical model calculations for stable isotope neutron capture*, Physical Review C **90**, 034619 (2014), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
* Mary Beard, Ethan Uberseder, and M. Wiescher, *Statistical Model Calculations for (n,g) Reactions*, EPJ Web of Conferences **93**, 04001 (2015), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
  
## Using Sapphire
If you are using *Sapphire* for your research, please send us a reference of your work and we will list it here.
* X. Fang, W. P. Tan, M. Beard, R. J. deBoer, G. Gilardy, H. Jung, Q. Liu, S. Lyons, D. Robertson, K. Setoodehnia, C. Seymour, E. Stech, B. Vande Kolk, M. Wiescher, R. T. deSouza, S. Hudan, V. Singh, X. D. Tang, and E. Uberseder, *Experimental measurement of 12C+16O fusion at stellar energies*, Physical Review C **96**, 045804 (2017), DOI: [10.1103/PhysRevC.96.045804](https://dx.doi.org/10.1103/PhysRevC.96.045804)
* M. Beard, E. Uberseder, R. Crowter, and M. Wiescher, *Comparison of statistical model calculations for stable isotope neutron capture*, Physical Review C **90**, 034619 (2014), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
* Mary Beard, Ethan Uberseder, and M. Wiescher, *Statistical Model Calculations for (n,g) Reactions*, EPJ Web of Conferences **93**, 04001 (2015), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
