# Atomistic kinetic Monte Carlo (A-kMC) model containting reactions that lead to the formation of defects from a pristine surface. 
### Aug 2020 - Dec 2021

The model is developed using separate, interacting lattices for the carbon and oxygen sites. The carbon lattice includes all possible carbon sites in a multilayer graphene structure representative of an HOPG material. The use of separate lattices makes it possible to effectively distinguish the oxygen and carbon sites and allow them to communicate throughout the simulation as various chemical interactions are tracked. 

## Required python packages
- numPy
- pandas
- cantera
- matplotlib
- sys

## Required data

Carbon_#x#_#.pkl: Dataframes containing necessary lists related to carbon lattice sites
Epoxies_#x#_#.pkl: Dataframes containing necessary lists related to epoxy lattice sites
coordinates_#_#.cfg: File containing indices and coordinates of carbon lattice sites 
coordinates_O_#_#.cfg: File containing indices and coordinates of epoxy lattice sites 

coordinates_#_#.cfg & coordinates_O_#_#.cfg for carbon and epoxy lattice were obtined as dump files from LAMMPS. 'in.carbon' & 'in.epoxy' files can generate these files for a given size and thickness of graphite. 

Python file KMC_initialization.py takes as input 'size' and 'no_of_layers' as input and outputs 'Carbon_#x#_#.pkl' and 'Epoxies_#x#_#.pkl' files (given that coordinates_#_#.cfg and coordinates_O_#_#.cfg are also present in the same folder).

## Description of scripts

### KMC.py
This script runs the A-kMC simulation with the following arguments:
_____________________________________
| Parameter | Description | Default |
|-----------|-------------|---------|
| `size` | size of graphene sheets | required |
| `no_of_layers` | no. of graphene layers | required |
| `temp` | temperature of gas [K] | required |
| `pressure` |pressure of gas mixture [Pa] | required |
| `temp_surface` | temperature of graphite surface [K] | `temp` |
| `save` | no. of iterations to save files | 1000000 |
| `walltime_max` | maximum walltime after which simulation stops [hrs] | 48 |
| `restart_iteration` | last run iteration to restart simulation from | -1 |

If a maximum walltime is specified, the script will run until it is reached and save a restart simulation file at the corresponding KMC iteration (to re-run the script from the last run KMC iteration). If a maximum walltime is not provided, no simulation file is saved at the corresponding iteration.

If a restart iteration is provided, KMC.py will run an A-kMC simulation starting at the state saved in './Dataframes/Df_restart_#Pa_#K_#_#_restart_iteration.pkl'

## How to use

- Clone this github repository
- Run KMC.py with the necessary arguments to run the A-kMC simulation.
- Run visualize.py with the necessary arguments to visualize simulation as LAMMPS dump files.
- Run statistics.py with the necessary arguments to get simulation statistics.

## Example

```bash
$ python KMC.py 20 5 2200 10000 --save 1000000 --walltime_max 4.5
```
This will generate 'Df_10000Pa_2200K_20_5_*.pkl' files in the Dataframes folder that will then be used for visualization and analysis.

To visualize these files run,

```bash
$ python visualize.py 20 5 2200 10000 --center 170
```
This will generate 'Sim_10000Pa_2200K_20_5_*pkl' files in the Results folder that can be visualized in Ovito.

To get simulation statistics run,

```bash
$ python visualize.py 20 5 2200 10000 > ./Results/output.txt
```
Details of the simulation such as adsorption and CO formaiton statistics will be recorded in the Results folder.

