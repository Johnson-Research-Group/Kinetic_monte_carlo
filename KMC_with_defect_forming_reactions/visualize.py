import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import glob
import re
import argparse

def parseArguments():
	
	"""
	This function creates argument parser
	"""

	parser = argparse.ArgumentParser()

	# required
	parser.add_argument("size", help="graphene sheet size", type=int)
	parser.add_argument("no_of_layers", help="no. of graphene sheets", type=int)
	parser.add_argument("temp", help="gas temperature", type=float)
	parser.add_argument("pressure", help="gas mixture pressure", type=float)

	parser.add_argument("-c", "--center", help="center sheets around this site", type=int, default=-1)

	# Parse arguments
	args = parser.parse_args()

	return args

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


# ==================================================================================================
#                                        Parameters
# ==================================================================================================

args = parseArguments()

size = args.size # Size unit of graphene sheets
no_of_layers = args.no_of_layers # No. of graphene layers
T_g = args.temp # Temperature of gas [K]
P_stag = args.pressure # Pressure of gas mixture [Pa]

# -------------------------------------- Print simulation details ---------------------------------

print("For size %d"%size)
print("For number of layers %d"%no_of_layers)
print("For temperature: %d"%T_g)
print("For Pressure: %d"%P_stag)

# ==================================================================================================
#                                       Reaction Categories
# ==================================================================================================

# Adsorption of oxygen on basal plane
Cat1=["Adsorption","Split O2 Adsorpiton"]

# Diffusion reactions
Cat2=["O2 Recombination & Desorption", "Epoxy diffusion Fast","Epoxy diffusion between layers Fast","Epoxy diffusion","Epoxy diffusion between layers", "1-Ep to 2-Ep","2-Ep to 3-Ep","1-Ep to 4-Ep", "1-Ep to Epoxy-Ether","Epoxy-Ether to 1-Ep","Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation)"]

# Surface reactions
Cat3=["Lactone-Ether CO formation","Ether-latone-ether CO formation","lactone-ether CO2 formation","ether-lactone-ether CO2 formation","lactone CO formation at VACANCY","lactone-epoxy CO formation at VACANCY","lactone CO2 formation at VACANCY","lactone-epoxy CO2 formation at VACANCY","Lactone-Ether to Ether-Lactone-Ether"]

# CO formation reactions at defect edge
Cat4_=["O2 ZigZag Adsorption","O2 Armchair Adsorption","O ZigZag Adsorption","O ArmChair Adsorption","O Carbonyl Adsorption","O2 Split Diffusion ZZ1","O2 Split Diffusion AC1", "O2 Split Diffusion AC2","O2 Split Diffusion ZZ2", "O2 desorption ZZ", "O2 desorption AC", "No Epoxy ZZ CO formation", "1 Epoxy ZZ CO formation", "2 Epoxy ZZ CO formation", "1 carbonyl AC CO formation", "2 carbonyl AC CO formation", "1 carbonyl Epoxy AC CO formation", "2 carbonyl Epoxy AC CO formation", "1 carbonyl DZZ CO formation", "2 carbonyl DZZ CO formation", "3 carbonyl DZZ CO formation", "SBC CO formation", "SBC Epoxy CO formation"] 

# Island removal reactions
add_Cat4=[]
for num in range(2,36):
	add_Cat4.append("%s Island Removal"%str(num)) 

Cat4= Cat4_+add_Cat4

# Combine categories into two categories, 'Category 1' & 'Category 2'
Category1 = Cat1+Cat2 # Adsorption & diffusion - deals with epoxy lattice sites
Category2 = Cat3+Cat4 # Surface reactions, edge adsorption and CO formation - deals with carbon lattice sites

# Combine categories into a single list 'Category'
Category = Category1+Category2

N1 = len(Category1)
N2 = len(Category2)

# ==================================================================================================
#                                           Load Lattices
# ==================================================================================================

# -------------------------------------- Load carbon lattice ---------------------------------------

f1=open("coordinates_%s_%d.cfg"%(str(size),no_of_layers),"r")
file=f1.readlines()

coord_x=[]
coord_y=[]
coord_z=[]

count=0
for line in file:
	count += 1
linesno=count

for line in file[9:count]:
	coord_x.append(float(line.split()[1]))
	coord_y.append(float(line.split()[2]))
	coord_z.append(float(line.split()[3]))

f1.close()

coordC_x = np.array(np.reshape(coord_x, (len(coord_x),1)))
coordC_y=np.array(np.reshape(coord_y, (len(coord_y),1)))
coordC_z=np.array(np.reshape(coord_z, (len(coord_z),1)))
layernno = np.array([no_of_layers-1-int(math.ceil(x/3.4)) for x in coord_z]) 
layernno = np.array(np.reshape(layernno, (len(layernno),1)))
id_no = np.array([i for i in range(0,len(coord_z))]) 
id_no = np.array(np.reshape(id_no, (len(id_no),1)))

# Coordinate array for carbon lattice (x coordinates, y coordinates, z coordinates, layer no., site id)
coordinates = np.copy(np.hstack((coordC_x, coordC_y, coordC_z, layernno, id_no)))
coordinates = np.stack( sorted(coordinates, key=lambda x: -x[2]), axis=0 )
for k in range(len(coordinates)):
	coordinates[k][4] = k

# Max length of epoxy lattice in x,y & z direction
max_C_x = max(coordC_x)[0]
max_C_y = max(coordC_y)[0]
max_C_z = max(coordC_z)[0]

# Surface area of graphene sheets in m^2
Surf_area = ((max_C_y*max_C_x)*(10**(-20)))

# ------------------------------------ Create carbonyl lattice ------------------------------------

# Replicate carbon lattice for carbonyl lattice (x coordinates, y coordinates, z coordinates, layer no.)
coordinates_O_2 = np.copy(np.hstack((coordC_x, coordC_y, coordC_z+1.275, layernno)))
coordinates_O_2 = np.stack( sorted(coordinates_O_2, key=lambda x: -x[2]), axis=0 )

# ---------------------------------- Create Cyclic-ether lattice ----------------------------------

# Replicate carbon lattice for cyclic-ether lattice 
coordinates_O_3 = np.copy(coordinates)

# -------------------------------------- Load epoxy lattice ---------------------------------------

f1=open("coordinates_O_%s_%d.cfg"%(str(size),no_of_layers),"r")
file=f1.readlines()

coord_x=[]
coord_y=[]
coord_z=[]

count=0
for line in file:
	count += 1
linesno=count

for line in file[9:count]:
	coord_x.append(float(line.split()[1]))
	coord_y.append(float(line.split()[2]))
	coord_z.append(float(line.split()[3]))

f1.close()

layernno_O = np.array([no_of_layers-1-int(math.ceil((x-1.29403)/3.4)) for x in coord_z]) 
layernno_O = np.array(np.reshape(layernno_O, (layernno_O.shape[0],1)))
coord_x=np.array(np.reshape(coord_x, (len(coord_x),1)))
coord_y=np.array(np.reshape(coord_y, (len(coord_y),1)))
coord_z=np.array(np.reshape(coord_z, (len(coord_z),1)))

# Coordinate array for epoxy lattice (x coordinates, y coordinates, z coordinates, layer no.)
coordinates_O_1 = np.copy(np.hstack((coord_x, coord_y, coord_z, layernno_O)))
coordinates_O_1 = np.stack(sorted(coordinates_O_1, key=lambda x: -x[2]), axis=0 )

# Max length of epoxy lattice in x & y direction
max_O_x = max(coord_x)[0]
max_O_y = max(coord_y)[0]

# Number of sites in Carbon, Epoxy, Carbonyl and Cyclic-ether lattice
sites = len(coordinates) #Carbon
sites_O_1 = len(coordinates_O_1) #Epoxy
sites_O_2 = len(coordinates_O_2) #Carbonyl
sites_O_3 = len(coordinates_O_3) #Cyclic-ether

# -------------------------------- centering poitns on one atom id --------------------------------

# indiex of site to center graphene sheet
center = args.center

if center!=-1:
	# To account for bonds on edge of graphene sheet
	max_C_x+=1.42
	max_C_y+=1.2297


	if coordinates[center[row_ind],0]<max_C_x/2:
		limit_x = coordinates[center[row_ind],0] + max_C_x/2
	else: 
		limit_x = coordinates[center[row_ind],0] - max_C_x/2

	if coordinates[center[row_ind],1]<max_C_y/2:
		limit_y = coordinates[center[row_ind],1] + max_C_y/2
	else: 
		limit_y = coordinates[center[row_ind],1] - max_C_y/2

	coordinates[:,0][coordinates[:,0]<limit_x]+= max_C_x
	coordinates_O_1[:,0][coordinates_O_1[:,0]<limit_x]+= max_C_x
	coordinates_O_2[:,0][coordinates_O_2[:,0]<limit_x]+= max_C_x
	coordinates_O_3[:,0][coordinates_O_3[:,0]<limit_x]+= max_C_x
	coordinates[:,1][coordinates[:,1]<limit_y]+= max_C_y
	coordinates_O_1[:,1][coordinates_O_1[:,1]<limit_y]+= max_C_y
	coordinates_O_2[:,1][coordinates_O_2[:,1]<limit_y]+= max_C_y
	coordinates_O_3[:,1][coordinates_O_3[:,1]<limit_y]+= max_C_y

# ==================================================================================================
#                                     Get simulation statistics
# ==================================================================================================


dfs = sorted(glob.glob("./Dataframes/Df_%dPa_%dK_%s_%d_*.pkl"%(P_stag,T_g,str(size), no_of_layers)), key=numericalSort)

for name in dfs:
	output = pd.read_pickle(name)

	# Load dataframe
	flag, flag_O_1,flag_O_2, flag_O_3,t, time_step, t_wall, changing_av_ad, changing_av_ep, site_O_before, site_before, site_ad, site_ep, edge_C, cov_C, cov_ep, counter_cat = ((output.tail(1)).values.tolist())[0]

	print("============================== KMC Iteration %d =================================="%time_step)

	# ---------------------------- Save as dump file for Ovito visualization --------------------------

	inds_to_keep = [i for i, x in enumerate(flag) if x !=['dg']]  
	sites_itr = len(inds_to_keep)

	inds_to_keep_O_1 = [i for i, x in enumerate(flag_O_1) if x !='na' and x !='pr' and x !='dg']   
	sites_itr_O_1 = len(inds_to_keep_O_1)

	inds_to_keep_O_2 = [i for i, x in enumerate(flag_O_2) if x !='na' and x !='pr' and x !='dg'] 
	sites_itr_O_2 = len(inds_to_keep_O_2)

	inds_to_keep_O_3 = [i for i, x in enumerate(flag_O_3) if x !='na' and x !='pr' and x !='dg'] 
	sites_itr_O_3 = len(inds_to_keep_O_3)

	colorC = [int(coordinates[i][3])+1 for i in inds_to_keep]

	color_O_1=[no_of_layers+1 if flag_O_1[i] == 'ep' else no_of_layers+2 if  flag_O_1[i] == 'la-eth' else no_of_layers+3 if  flag_O_1[i] =='eth-la-eth' else no_of_layers+6 if  flag_O_1[i] =='ep-eth' else no_of_layers+7 if flag_O_1[i] == 'la-vac' else no_of_layers+9 for i in inds_to_keep_O_1]
	color_O_2=[no_of_layers+2 if flag_O_2[i] == 'la-eth' else no_of_layers+3 if  flag_O_2[i] =='eth-la-eth' else no_of_layers+4 if flag_O_2[i] == 'carb_O' else no_of_layers+7 if flag_O_2[i] == 'la-vac' else no_of_layers+8 if flag_O_2[i] == 'carb_O2' else no_of_layers+9 for i in inds_to_keep_O_2]
	color_O_3=[no_of_layers+5 if flag_O_3[i] == 'cyc-eth' else no_of_layers+9 for i in inds_to_keep_O_3]

	f = open("./Results/Sim_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step), "w")
	f.write("ITEM: TIMESTEP\n")
	f.write("%d\n"%offset_step)
	f.write("ITEM: NUMBER OF ATOMS\n")
	f.write("%d\n"%(sites_itr+sites_itr_O_1+sites_itr_O_2+sites_itr_O_3))#+sites_itr_O))
	f.write("ITEM: BOX BOUNDS pp pp ff\n")
	f.write("-1 171\n")
	f.write("-1 150\n")
	f.write("-1 100\n")
	f.write("ITEM: ATOMS id type x y z\n")
	for i in inds_to_keep:
		f.write("%d %d %.9f %.9f %.9f\n"%(i,colorC[inds_to_keep.index(i)],coordinates[i][0],coordinates[i][1],coordinates[i][2]))
	for i in inds_to_keep_O_1:
		f.write("%d %d %.9f %.9f %.9f\n"%(i+sites,color_O_1[inds_to_keep_O_1.index(i)],coordinates_O_1[i][0],coordinates_O_1[i][1],coordinates_O_1[i][2]))
	for i in inds_to_keep_O_2:
		f.write("%d %d %.9f %.9f %.9f\n"%(i+sites+sites_O_1,color_O_2[inds_to_keep_O_2.index(i)],coordinates_O_2[i][0],coordinates_O_2[i][1],coordinates_O_2[i][2]))
	for i in inds_to_keep_O_3:
		f.write("%d %d %.9f %.9f %.9f\n"%(i+sites+sites_O_1+sites_O_2,color_O_3[inds_to_keep_O_3.index(i)] ,coordinates_O_3[i][0],coordinates_O_3[i][1],coordinates_O_3[i][2]))
	f.close()
