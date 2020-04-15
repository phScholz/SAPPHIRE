# Sapphire - a statistical nuclear reaction and decay code

- [General Information](#general-information)
  - [Sapphire](#sapphire)
  - [Former Contributors](#former-contributors)
  - [License](#license)
- [Compilation & Execution](#compilation--execution)
  - [Requirements](#requirements)
  - [Build](#build)
  - [Run](#run)
- [Modules](#modules)
    - [reaction](#reaction)
    - [decayer](#decayer)
    - [old](#old)
    - [help](#help)
- [References](#references)
  - [About Sapphire](#about-sapphire)
  - [Using Sapphire](#using-sapphire)

--- 

## General Information

### Sapphire

* is an independent code which can calculate both radiative capture and inverse cross sections.

* is a new code, written in C++, designed primarily to simulate the decay of excited nuclei via the Monte Carlo technique. 

* generates realistic particle and γ-ray distributions resulting from the statistical decay of a compound nucleus having a predefined excitation energy,with decay probabilities based on the most recent experimental data compilations and state-ofthe-art theoretical models. The code is designed for large-scale multiprocessing and therefore is ideal for the generation of the large theoretical data sets often needed to determine experimental response functions without the need of excessive approximations. The transmission functions needed to form the probability distribution functions forMonte Carlo are identical to those required to calculate HF cross sections; therefore, the code can be used to calculate astrophysical reaction rates. 

* considers a range of projectile nuclei, and include four exit channels, namely: (particle,n), (particle,p), (particle,α), and (particle,γ ). 

* has the additional ability to calculate cross sections for photon induced reactions. 

* when used for cross-section calculations, does not include preequilibrium reactions, which makes them well suited to astrophysical studies because in these environments incident particles have energies that are typically less than approximately 10 MeV.

### Former Contributors
* original version of Sapphire by [Mary Beard](https://isnap.nd.edu/people/group-pages/mary-beard/). 
* Ethan Uberseder

### License
Sapphire is currently published under the [GNU General Public License](GPL3.md).

## Compilation & Execution
### Requirements
Compiling Sapphire requires 
* [CMake](https://cmake.org/),
* [ROOT](https://root.cern.ch/),
* [boost](https://www.boost.org/) and
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). 

Check the find_packages in `cmake/Modules` whether cmake is looking for these packages in the right place if cmake cannot find them.

For an `automized generation of documentation files`, one needs [doxygen](http://doxygen.nl/) installed on the system. If its not installed, documentation is not automatically generated. 

### Build
#### Automatically
A simple bash script has been written, does the steps below automatically for you. Since it uses relative paths you need to be in your Sapphire main directory and then simply execute:
```bash
bash build.sh
```

To delete the directories which have been created during the build process, simply run the following command in your Sapphire main directory:
```bash
bash clean.sh
```

#### Manually
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

### Run

To check if the code is running, simply type 

```
./sapphire old X+a
```

or 

```
./sapphire reaction X+a
```

where X is the target nucleus and a is the projectile. Examples are 25Mg+a or 60Fe+n.  

---

## Modules

The new version of the Sapphire code is modularized to increase the maintainability and readability of the code basis. At the moment there the modules [reaction](#reaction), [decayer](#decayer), and [old](#old) available. See below for a more detailed explanation.

### reaction
The reaction module is responsible for the calculation of reaction cross sections and reaction rates. 

For quick calculations, one can simply type

```
./sapphire reaction X+a
```

where X+a is a reaction string with X the target nucleus and a the projectile. This command starts a default calculation, i.e. all input parameters are set to default.

For a more detailed calculations, one should use the possibility to parse an input file to sapphire by typing 

```
./sapphire reaction inputFile.ini
```

The file `inputFile.ini` will be parsed by the [Boost.PropertyTree](https://www.boost.org/doc/libs/1_72_0/doc/html/property_tree.html) library. The gneral format and keyword descriptions can be found [here](INPUTFILE.md).


### decayer
Not yet fully implemented.

### old
This module provides the same functionallity as the original Sapphire code. Everything stays the same except the keyword ```old```. 
Thus instead of running ```./sapphire 60Fe+a``` one needs to run ```./sapphire old 60Fe+a```.

In future more and more of the features of the original Sapphire code will be reimplemented as modules. This part of the code can still be used but won't be maintained in future.

### help
Print the help text for the different modules.

---

## References
### About Sapphire
If you are using *Sapphire* you should cite one or more of the references below.

* M. Beard, E. Uberseder, R. Crowter, and M. Wiescher, *Comparison of statistical model calculations for stable isotope neutron capture*, Physical Review C **90**, 034619 (2014), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
* Mary Beard, Ethan Uberseder, and M. Wiescher, *Statistical Model Calculations for (n,g) Reactions*, EPJ Web of Conferences **93**, 04001 (2015), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
  
### Using Sapphire
If you are using *Sapphire* for your research, please send us a reference of your work and we will list it here.
* X. Fang, W. P. Tan, M. Beard, R. J. deBoer, G. Gilardy, H. Jung, Q. Liu, S. Lyons, D. Robertson, K. Setoodehnia, C. Seymour, E. Stech, B. Vande Kolk, M. Wiescher, R. T. deSouza, S. Hudan, V. Singh, X. D. Tang, and E. Uberseder, *Experimental measurement of 12C+16O fusion at stellar energies*, Physical Review C **96**, 045804 (2017), DOI: [10.1103/PhysRevC.96.045804](https://dx.doi.org/10.1103/PhysRevC.96.045804)
* M. Beard, E. Uberseder, R. Crowter, and M. Wiescher, *Comparison of statistical model calculations for stable isotope neutron capture*, Physical Review C **90**, 034619 (2014), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
* Mary Beard, Ethan Uberseder, and M. Wiescher, *Statistical Model Calculations for (n,g) Reactions*, EPJ Web of Conferences **93**, 04001 (2015), DOI: [10.1103/PhysRevC.90.034619](https://dx.doi.org/10.1103/PhysRevC.90.034619)
