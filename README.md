# Hyperthermal Oxidation of Graphite using Atomistic kinetic Monte Carlo (A-kMC) Simulations
An atomistic kinetic Monte Carlo model is used to simulate defect generation in multi-layer graphene due to atomic and molecular oxygen. The purpose of these simulations is to study the erosion of Thermal Protection System (TPS) materials at high temperatures during space reentry to earth. 

Simulations based on this KMC model reveal details of the pitting process driven by oxygen adsorption and surface diffusion in multi-layer graphene for specific temperatures and partial pressures of atomic and molecular oxygen. 

The model implemented here consists of 44 reactions between carbon and oxygen, categorized as adsorption, diffusion, surface and gasification reactions. The following figures shows the oxidation of graphite consisting of five graphene layers at 1700 K and 10kPa.

<div style="text-align: center;">
<img src="Images/oxidation_1700K_10kPa.png" width="600" height="300">
</div>

As the temperature increases from 1300 K to 2200 K, the shape of the defects transitions from hexagonal to branched. This change is attributed to the preferential desorption of carbonyl groups from certain edge sites as the temperature rises.

The following figures show the change in defect shape form 1300 K to 2200 K
![Image](Images/temperatures.jpg)

The A-kMC model is useful for investigating reaction kinetics and its behavior on defect formation at the atomic scale.

Please refer to the citations listed for more information.

## Projects in this repository
- Atomistic-kinetic Monte Carlo (A-kMC), with defect forming reactions
- Atomistic-kinetic Monte Carlo (A-kMC), without defect forming reactions
- Direct Simulation Monte Carlo coupled with Atomistic-kinetic Monte Carlo

### References
[1] Edward, Sharon, and Harley Johnson. "Atomistic Multi-Lattice Kinetic Monte Carlo (KMC) Modeling of Hyperthermal Oxidation of Multi-Layer Graphene." AIAA SCITECH 2022 Forum. 2022.

[2] Edward, Sharon, and Harley T. Johnson. "Oxidation of Multi-layer Graphite Using an Atomistic Multi-lattice Kinetic Monte Carlo Model." The Journal of Physical Chemistry C 127.34 (2023): 16938-16949.

[3] Edward, Sharon, Moon-ki Choi, and Harley T. Johnson. "Application of Time Series Analysis on Kinetic Monte Carlo Simulations of Hyperthermal Oxidation of Graphite." The Journal of Physical Chemistry C (2025).