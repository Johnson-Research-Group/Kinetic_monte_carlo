## Atomistic kinetic Monte Carlo (A-kMC) model containting reactions that lead to the formation of defects from a pristine surface. 
# Aug 2020 - Dec 2021

The model is developed using separate, interacting lattices for the carbon and oxygen sites. The carbon lattice includes all possible carbon sites in a multilayer graphene structure representative of an HOPG material. The use of separate lattices makes it possible to effectively distinguish the oxygen and carbon sites and allow them to communicate throughout the simulation as various chemical interactions are tracked. 

First run KMC_initialization.py to get the 'Carbon_#x#_#.pkl' and 'Epoxies_#x#_#.pkl'

Make sure

Next run KMC.py:

This script is set to run at:
Gas and surface temperature: 1500K (Line 834-835)
Stagnation pressure: 10000 Pa (Line 850)
Mole fraciton of atomic oxygen: 5.68167898e-06 (Obtained from NASA CEA database)
Mole Fraction of molecular oxygen: 1.93874074e-01 (Obtained from NASA CEA database)
Partial pressure of atomic oxygen = Stagnation pressure * Mole fraciton of atomic oxygen (Line 851)
Partial pressure of molecular oxygen = Stagnation pressure * Mole fraciton of molecular oxygen (Line 852)

Before running the script, set initial parameters in lines 653-654
Line 653: size of graphene sheets (as set in KMC_initialization.py)
Line 737: number of graphene layers

Next run write_to_files.py:

Before running the script, set initial parameters in lines 13-16
Line 13: [List of all sizes you want to run]
Line 14: Number of layers of graphene for all sizes
Line 15: [List of final timestep saved for respective size]
Line 16: [List of total number of carbon atoms for respective size]

The coordinates_#_#.cfg, coordinates_O_#_#.cfg, Carbon_#x#_#.pkl and Epoxies_#x#_#.pkl files need to be in the same folder while running the scripts.
coordinates_#_#.cfg and coordinates_O_#_#.cfg files were obtained from LAMMPS to load coordinate points for carbon and oxygen lattice points respectively.
