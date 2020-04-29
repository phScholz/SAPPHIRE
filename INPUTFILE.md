# Inputfile description and keywords - contents

1. [Keywords Table](#keywords-table)
2. [Keywords](#keywords)
3. [Example Inputfile](#example-inputfile)

---

# Keywords Table
## General

| Keyword                       | Short explanation                                                      |
| :---------------------------- | :--------------------------------------------------------------------- |
| [MassTable](#masstable)       | Path to the mass tables                                                |
| [GDRParams](#gdrparams)       | Path to the file with gdr parameters                                   |
| [LevelDir](#leveldir)         | Path to the directory which contains the level files                   |
| [SpinFile](#spinfile)         | Path to the file which contains the spin information                   |
| [ProtonModel](#threads)       | Sets the model for proton transmission function.                       |
| [NeutronModel](#threads)      | Sets the model for neutron transmission function.                      |
| [GammaModel](#threads)        | Sets the model for gamma transmission function.                        |
| [AlphaModel](#threads)        | Sets the model for alpha transmission function.                        |
| [LevelDensity](#LevelDensity) | Sets the model for the Nuclear Level Density.                          |
| [cTable](#cTable)             | Value to normalize the level spacing at the neutron separation energy. |

## CrossSection
| Keyword              | Short explanation                                                         |
| :------------------- | :------------------------------------------------------------------------ |
| Reaction             | Reaction string,e.g., 25Mg+a                                              |
| EnergyFile           | Path to the file containing a list of energies                            |
| ReactionFile         | Path to the file containing a list of reactions                           |
| CalcRates            | Turn calculation of rates on or off (if neutron)                          |
| CalcAverageWidth     | Calculates the average s-wave radiative width at threshold and exits.     |
| ResidualGamma        | Toggles if residual or total capture cross section is calculated.         |
| ResidualNeutron      | Toggles if neutron residual or total capture cross section is calculated. |
| ResidualProton       | Toggles if proton residual or total capture cross section is calculated.  |
| ResidualAlpha        | Toggles if alpha residual or total capture cross section is calculated.   |
| CalculateGammaCutoff | Toggles if Gamma-Cutoff energy should be calculated or not.               |
| EntranceState        | Give the number of the initial state, e.g. 0 for groundstate              |
| g_ExitStates         | Define how many partial cross section should be calculated.               |
| n_ExitStates         | Define how many partial cross section should be calculated.               |
| p_ExitStates         | Define how many partial cross section should be calculated.               |
| a_ExitStates         | Define how many partial cross section should be calculated.               |

## Decay
| Keyword           | Short explanation                              |
| :---------------- | :--------------------------------------------- |
| [Suffix](#suffix) | Suffix which can be given to the output files. |
| PreEqConf         | Sets the initial exciton configuration.        |

## Random


# Keywords
## General
### MassTable
The default mass tables can be found in `/tables/masses/`. Currently, Sapphire uses the mass excesses from the [Reference Input Parameter Library (RIPL-3)](https://www-nds.iaea.org/RIPL-3/) to calculate atomic masses and Q-Values. The experimental mass excesses were evaluated by [G. Audi and A.H. Wapstra](https://www.sciencedirect.com/science/article/abs/pii/0375947495004459). 

If no experimental mass excesses are available, Sapphire uses the ground state properties of the [Finite Range Droplet Model (FRDM)](https://www.sciencedirect.com/science/article/pii/S0092640X85710029) from RIPL-3.

Users who want to use different ground state properties are free to edit the data tables or exchange them with other ones. However, the file structure needs to stay the same as defined by RIPL-3.

### GDRParams

### LevelDensity
Currently available models:

0. **Rauscher Level Density**
1. **Microscopical Level Density HFB-BSk14** 
   
   Quotation from RIPL-3: 

   "The nuclear level density is coherently obtained on the basis of the single-particle level scheme and pairing energy derived at the ground state deformation based on the BSk14 Skyrme force ([S. Goriely, M. Samyn, J.M. Pearson, Phys. Rev. C75 (2007) 064312](https://doi.org/10.1103/PhysRevC.78.064307)). 
   
   Additionally, the phenomenological level density parameters ctable and ptable are tabulated in files (*.cor) by fitting the HFB calculated curve to the RIPL II recommended spacings of s-wave neutron resonances D0 and to the cumulative number of low-lying levels ([S. Goriely, S. Hilaire, A.J. Koning, Phys. Rev. C *78* (2008) 064307](http://dx.doi.org/10.1103/PhysRevC.78.064307))".



### cTable

Normalization factor $c$ for level densities which are based on LevelDensityTables:
$$\rho(E,J,P)_{norm} = e^{c\sqrt{E}}\rho(E,J,P)$$

With this parameter the level density can be fitted to level spacing at the neutron separation energy $S_n$.


### LevelDir
The known levels are also derived from RIPL-3 and can be found in `/tables/levels/`. 
### SpinFile
### Suffix

## CrossSection

## Decay


# Example Inputfile
Below an example of an input file for the reaction module is shown.
The paths to the input files for masses, levels, etc. are created automatically from the source directory during the building process of the code. These only have to been edited, if the source directory has been moved since the last compilation.


```
[General]
MassTable            = ./SAPPHIRE/tables/masses.dat
GDRParams            = ./SAPPHIRE/tables/ripl3_gdr_parameters.dat
Leveldir             = ./SAPPHIRE/levels/
SpinFile             = ./SAPPHIRE/tables/spinod.dat
Suffix               = 

[CrossSection]
Reaction             = 25Mg+a
EnergyFile           = 
ReactionFile         = ./SAPPHIRE/examples/reactionFile
CalcRates            = 0
CalcAverageWidth     = 0
ResidualGamma        = 1
ResidualNeutron      = 0
ResidualProton       = 0
ResidualAlpha        = 0
CalculateGammaCutoff = 0
EntranceState        = 0
g_ExitStates         = -1
n_ExitStates         = -1
p_ExitStates         = -1
a_ExitStates         = -1
```