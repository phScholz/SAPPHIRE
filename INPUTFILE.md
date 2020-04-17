# Inputfile description and keywords - contents

- [Keywords Table](#keywords-table)
- [Keywords](#keywords)
- [Example Inputfile](#example-inputfile)

---

# Keywords Table
## General

| Keyword       | Short explanation 
| :--- | :---      
| [MassTable](#masstable) | Path to the mass tables   
| [GDRParams](#gdrparams) | Path to the file with gdr parameters
| [LevelDir](#leveldir)   | Path to the directory which contains the level files
| [SpinFile](#spinfile)   | Path to the file which contains the spin information
| [ProtonModel](#threads)  | Sets the model for proton transmission function.
| [NeutronModel](#threads) | Sets the model for neutron transmission function.
| [GammaModel](#threads)   | Sets the model for gamma transmission function.
| [AlphaModel](#threads)   | Sets the model for alpha transmission function.

## CrossSection
| Keyword       | Short explanation 
| :--- | :---      
|Reaction            | Reaction string,e.g., 25Mg+a 
|EnergyFile          | Path to the file containing a list of energies
|ReactionFile        | Path to the file containing a list of reactions
|CalcRates           | Turn calculation of rates on or off (if neutron)
|CalcAverageWidth    | Calculates the average s-wave radiative width at threshold and exits.
|ResidualGamma       | Toggles if residual or total capture cross section is calculated.
|ResidualNeutron     | Toggles if neutron residual or total capture cross section is calculated.
|ResidualProton      | Toggles if proton residual or total capture cross section is calculated.
|ResidualAlpha       | Toggles if alpha residual or total capture cross section is calculated.
|CalculateGammaCutoff| Toggles if Gamma-Cutoff energy should be calculated or not.
|EntranceState       | Give the number of the initial state, e.g. 0 for groundstate
|g_ExitStates        | Define how many partial cross section should be calculated.
|n_ExitStates        | Define how many partial cross section should be calculated.
|p_ExitStates        | Define how many partial cross section should be calculated.
|a_ExitStates        | Define how many partial cross section should be calculated.

## Decay
| Keyword       | Short explanation 
| :--- | :---      
| [Suffix](#suffix)       | Suffix which can be given to the output files.
| PreEqConf | Sets the initial exciton configuration.

# Keywords
## General
### MassTable
### GDRParams
### LevelDir
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
