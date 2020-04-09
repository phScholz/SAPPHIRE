This is the original version of Sapphire by [Mary Beard](https://isnap.nd.edu/people/group-pages/mary-beard/). 
Only the README-file has been changed.


Compiling Sapphire requires CMake,ROOT,and GSL.  

To compile:

1. Create a build directory under the main Sapphire directory 
  (i.e. mkdir build), and change to that directory (i.e. cd build/).
2. Run CMake against the main directory, optionally specifying the 
  desired C++ compiler (i.e. cmake -DCMAKE_CXX_COMPILER=icpc  -DCMAKE_C_COMPILER=icc ..).
3. Type make install to build Sapphire.  The executable is put in
  the build directory.

Be aware that Sapphire links the paths to the needed tables at compile time.
While the executable can be moved, the main Sapphire directory should stay in
place.  If moved, the build process should be repeated.

To execute the code just enter

sapphire X+a

where X is the heavy nucleus and a is the projectile. Examples are 25Mg+a or 60Fe+n.  
