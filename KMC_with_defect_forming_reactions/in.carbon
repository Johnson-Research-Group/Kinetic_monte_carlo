
#-------------Units, Dimensions------------
units	 real

dimension 3

newton on

processors * * *

boundary p p f

atom_style charge

atom_modify map yes

#----------Create atoms---------------

variable l equal 1.42
variable s equal 3.4
variable a1 equal 3*$l
variable a2 equal (sqrt(3))*$l
variable a3 equal $s*2

variable a equal 1/6
variable b equal 2/3
variable c equal 1/3
variable d equal 5/6

# ====================== PARAMETERS FOR GRAPHiTE ========================
variable u equal 20 # graphene sheet size
variable lay equal 15 # thickness of epoxy lattice in Angstroms (layers in epoxy lattice should match that of carbon lattice)

lattice custom 1.0 a1 ${a1} 0.0 0.0 a2 0.0 ${a2} 0.0 a3 0.0 0.0 ${a3} &
                    basis 0 0 0 basis $a 0.5 0 basis 0.5 0.5 0 basis $b 0 0 basis 0 0 0.5 basis $c 0 0.5 basis 0.5 0.5 0.5 basis $d 0.5 0.5

variable b1 equal ${a1}*${u}*2-0.71
variable b2 equal ${a2}*${u}*3
variable b3 equal ${a3}
variable b4 equal ${b3}*2
variable ho equal 2*${b3}+15

region mybox block -0.72 ${b1} 0.0 ${b2} -35.0 100.0 units box
create_box 2 mybox 
region C_domain_1 block -0.72 ${b1} 0.0 ${b2} 0 14 units box #8 for 3 layers 14 for 5 3 for 1 layer 4 for 2 layers
create_atoms 1 region C_domain_1

group carbon type 1
group oxygen type 2

dump graph all custom 25 coordinates.cfg id x y z 
#--------------Mass and Velocity--------

mass 1 12.0107
mass 2 15.999

#---------Neighbor----------------------

neighbor    2 bin
neigh_modify    every 1 delay 0 check yes

#------------Timestep--------------------

timestep 0.25 #femtoseconds

#---------Pair Style & Coeff------------

pair_style reax/c NULL 
pair_coeff * * ffield.reax.cho C O

#--------------Fix Commands-------------

fix qeq_reax all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#------------Energy minimization-----

fix relax carbon box/relax x 0.0 y 0.0 couple none fixedpoint 0 0 0

compute ctemp carbon temp
compute peratom carbon pe/atom
compute pe carbon reduce sum c_peratom
compute otemp oxygen temp
compute temp_com carbon temp/com # Temperature of all Carbon layers without center-of-mass velocity

thermo_style custom step etotal c_pe c_otemp c_ctemp lx ly
thermo_modify lost ignore
thermo 100

min_style cg
min_modify line quadratic
minimize 1.0e-30 1.0e-30 100000 1000000

fix dynamics all nve

run 0

