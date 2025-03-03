# Atomistic kinetic Monte Carlo (A-kMC) model containting reactions that lead to the formation of defects from a pristine surface. 
### Aug 2020 - Dec 2021

The model is developed using separate, interacting lattices for the carbon and oxygen sites. The carbon lattice includes all possible carbon sites in a multilayer graphene structure representative of an HOPG material. The use of separate lattices makes it possible to effectively distinguish the oxygen and carbon sites and allow them to communicate throughout the simulation as various chemical interactions are tracked. 

## Required python packages
- NumPy
- Pandas
- Cantera
- Matplotlib
- Sys

## Required data

Carbon_#x#_#.pkl: Dataframes containing necessary lists related to carbon lattice sites
Epoxies_#x#_#.pkl: Dataframes containing necessary lists related to epoxy lattice sites
coordinates_#_#.cfg: File containing indices and coordinates of carbon lattice sites 
coordinates_O_#_#.cfg: File containing indices and coordinates of epoxy lattice sites 

coordinates_#_#.cfg & coordinates_O_#_#.cfg for carbon and epoxy lattice were obtined as dump files from LAMMPS. 'in.carbon' & 'in.epoxy' files can generate these files for a given size and thickness of graphite. 

Python file KMC_initialization.py takes as input 'size' and 'no_of_layers' as input and outputs 'Carbon_#x#_#.pkl' and 'Epoxies_#x#_#.pkl' files (given that coordinates_#_#.cfg and coordinates_O_#_#.cfg are also present in the same folder).

## How to use

- Ensure that 'Carbon_#x#_#.pkl', 'Epoxies_#x#_#.pkl', 'coordinates_#_#.cfg' and 'coordinates_O_#_#.cfg' are in the same folder with pyhton scripts.
- Run KMC.py '# graphene sheet size (int)' '# no. of graphene layers (int)' 'maximum iterations (int)' 'iterations to save files (int)' 'maximum walltime (int)'
- Run write_to_files.py to visualize simulation as LAMMPS dump files

## Example

```bash
# This is a terminal command
$ echo "Hello, World!"
Hello, World!
