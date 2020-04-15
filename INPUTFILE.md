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
| [Suffix](#suffix)       | Suffix which can be given to the output files

## CrossSection
| Keyword       | Short explanation 
| :--- | :---      
|Reaction            | 
|EnergyFile          | 
|ReactionFile        | 
|CalcRates           | 
|CalcAverageWidth    | 
|ResidualGamma       | 
|ResidualNeutron     | 
|ResidualProton      | 
|ResidualAlpha       | 
|CalculateGammaCutoff| 
|EntranceState       | 
|g_ExitStates        | 
|n_ExitStates        | 
|p_ExitStates        | 
|a_ExitStates        | 

## Decay

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