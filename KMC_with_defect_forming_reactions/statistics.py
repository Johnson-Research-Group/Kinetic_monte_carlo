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

	# Parse arguments
	args = parser.parse_args()

	return args

def numericalSort(value):

	
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def what_edge_carbon_is_it(site):
	list_n =[]
	if edge_C[site] =='edg':
		list_n = [l for l in near_neigh[site] if flag[l] != ['dg']]
		if len(list_n) == 2:
			if [flag[p] for p in near_neigh[list_n[0]]].count(['dg']) == [flag[p] for p in near_neigh[list_n[1]]].count(['dg']) == 0 :
				return ["ZZ",list_n]
			elif [flag[p] for p in near_neigh[list_n[0]]].count(['dg']) == [flag[p] for p in near_neigh[list_n[1]]].count(['dg']) == 1:
				return ["DZZ",list_n]
			elif [flag[p] for p in near_neigh[list_n[0]]].count(['dg']) == [flag[p] for p in near_neigh[list_n[1]]].count(['dg']) == 1:
				return ["TZZ",list_n]
			elif [flag[p] for p in near_neigh[list_n[0]]].count(['dg']) != [flag[p] for p in near_neigh[list_n[1]]].count(['dg']) and any([flag[p] for p in near_neigh[i]].count(['dg']) == 1 for i in list_n)==True:
				return ["AC",list_n]
			else:
				return ["ZZ",list_n]
		elif len(list_n) == 1:
			return ["SBC",list_n] #SBC - Single bonded carbon
		else:
			return ["Island",list_n]
	else:
		return ["Not in edge",list_n]

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
print()

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

# ==================================================================================================
#                                    Load Carbon & Epoxy lists
# ==================================================================================================

# Load dataframes
output1 = pd.read_pickle("Carbon_%dx%d_%d.pkl"%(size,size,no_of_layers))
output2 = pd.read_pickle("Epoxies_%dx%d_%d.pkl"%(size,size,no_of_layers))

# Assign loaded lists for carbon sites
near_neigh = list(output1['Near Neighbors'])
epoxy_list = list(output1['Epoxy List'])
check_above = list(output1['Check Above'])
check_below = list(output1['Check Below'])
ZigZag_neigh = list(output1['ZigZiag Neighbors'])
Arm_chair_neigh = list(output1['ArmChair Neighbors'])
near_neigh_C_O_7 = list(output1['Near Neighbors C O 7'])
near_neigh_C_8 = list(output1['Near Neighbors C 8'])

# Assign loaded lists for epoxy sites
epoxy_O = list(output2['Epoxy O'])
epoxy_C_move = list(output2['Epoxy C Move'])
epoxy_O_move = list(output2['Epoxy O Move'])
epoxy_O_orient = list(output2['Epoxy O Orient'])
epoxy_O_split = list(output2['Epoxy O Split'])
near_neigh_8 = list(output2['Near Neighbors 8'])
near_neigh_O_7 = list(output2['Near Neighbors O 7'])

# ==================================================================================================
#                                     Record simulation statistics
# ==================================================================================================

physical_time =[]
carbon_removed = []
basal_oxygen = []
defect_edge_oxygen = []

dfs = sorted(glob.glob("./Dataframes/Df_%dPa_%dK_%s_%d_*.pkl"%(P_stag,T_g,str(size), no_of_layers)), key=numericalSort)

for name in dfs:
	output = pd.read_pickle(name)

	flag, flag_O_1, flag_O_2, flag_O_3, t, time_step, t_wall, edge_C, cov_C,cov_ep,counter_cat = ((output.tail(1)).values.tolist())[0]

	# ----------------------------------------- Print times ------------------------------------------

	print("---------------------------------------------------------------------------------")
	print("                              KMC Iteration %d                                 "%time_step)
	print("---------------------------------------------------------------------------------")

	print()
	print("Physical Time: %s secs"%"{:.2e}".format(t))
	print("Wall Time: %.7f secs"%(t_wall))
	print()
	print()
	# ----------------------------------- Adsorbed oxygen statistics ----------------------------------

	total_ads_O = counter_cat[0]+counter_cat[24]+counter_cat[25]+counter_cat[26]
	total_ads_O2 = counter_cat[1]+counter_cat[22]+counter_cat[23]
	total_ads = total_ads_O + 2*total_ads_O2

	total_basal_des = counter_cat[2]
	total_edge_O2_des = counter_cat[31] + counter_cat[32]
	total_des = total_basal_des + 2*total_edge_O2_des

	total_ads_basal = counter_cat[0]+2*counter_cat[1] - total_basal_des
	total_ads_edge = counter_cat[24]+counter_cat[25]+counter_cat[26]+2*(counter_cat[22]+counter_cat[23]-total_edge_O2_des)

	print(" ==================== Adsorption & desorption statistics ========================= ")

	print()
	print("Total number of oxygen adsorbed: %d"%total_ads)
	print("Total number of atomic oxygen adsorbed: %d"%total_ads_O)
	print("Total number of molecular oxygen adsorbed: %d"%total_ads_O2)
	print()
	print("Total number of oxygen desorbed: %d"%total_des)
	print("Total number of epoxies desorbed: %d"%total_basal_des)
	print("Total number of molecular oxygen desorbed from edge: %d"%total_edge_O2_des)
	print()
	print("Total effectively adsorbed on the basal plane (subtracting desorbed oxygen): %d"%total_ads_basal)
	print("Total effectively adsorbed on defect edges (subtracting desorbed oxygen): %d"%total_ads_edge)
	print()
	print()
	basal_oxygen.append(total_ads_basal)
	defect_edge_oxygen.append(total_ads_edge)

	# -------------------------------------- CO formation statistics -----------------------------------

	print(" ======================== CO formation statistics ================================ ")

	print()
	total_removed = flag.count('dg')
	total_removed_CO2 = sum([counter_cat[i] for i in [15,16,19,20]])
	total_removed_CO = sum([counter_cat[i] for i in [13,14,17,18]+list(range(33,45))])
	print("Total number of carbon atoms removed: %d"%total_removed)
	print("Total CO2 formation: %d"%total_removed_CO2)
	print("Total CO formation: %d"%total_removed_CO)

	ZZ_CO =0
	AC_CO =0
	DZZ_CO =0
	SBCC_CO =0

	for i in [21,22,23]:
		i=i+12
		ZZ_CO += counter_cat[i]

	print("Number of zigzag sites removed: %d"%ZZ_CO)

	for i in [24,25,26,27]:
		i=i+12
		AC_CO += counter_cat[i]

	print("Number of armchair sites removed: %d"%AC_CO)

	for i in [28,29,30]:
		i=i+12
		DZZ_CO += counter_cat[i]

	print("Number of double zigzag sites removed: %d"%DZZ_CO)

	for i in [31, 32]:
		i=i+12
		SBCC_CO += counter_cat[i]

	print("Number of single bonded carbon sites removed: %d"%SBCC_CO)

	print("Total number of carbon atoms removed as islands: %d"%(total_removed-total_removed_CO2-total_removed_CO))
	print()
	print()

	# ----------------------------------- Print number of edge sites ----------------------------------

	print(" ============================= Other statistics ================================== ")

	print()
	print("Number of carbon edge sites: %d"%edge_C.count('edg'))
	print("Number of exposed carbon sites: %d"%cov_C.count('uncov'))
	edge_type = []
	for i in range(sites):
		if edge_C[i]=='edg':
			typ = what_edge_carbon_is_it(i)
			if typ[0]=="ZZ":
				edge_type.append('z')
			elif typ[0]=="AC":
				edge_type.append('a')
			else:
				pass
	zigzag = edge_type.count('z')
	armchair = edge_type.count('a')

	if armchair!=0:
		print("Ratio of zigzag to armchair sites: %.4f"%(zigzag/armchair))
	else:
		print("No armchair sites, number of zigzag sites: %d"%zigzag)
	print("Number of initiated defects: %d"%sum([counter_cat[i] for i in range(13,21)]))
	print()
	print()

	carbon_removed.append(total_removed)
	physical_time.append(t)

# -------------------------- Print number of occurences for each reaction -------------------------

print(" ============ Number of occurences for each reaction ============ ")

print()
print("No.| Category | Number of occurences")
for cat in range(len(Category)):
	print(cat,"|", Category[cat],"|", counter_cat[cat])
print()
print()

# -------------------------------------- Rate of CO formation -------------------------------------

length = len(carbon_removed)-1
rate = carbon_removed[length]/physical_time[length]
print("Rate of carbon removal: %s [C atoms/s]"%"{:.2e}".format(rate))
print()

# ------------------------------------------- Save plot -------------------------------------------

# Convert time to milliseconds
physical_time = [i*10**3 for i in physical_time]

# Save plot
plt.figure()
plt.plot(carbon_removed, physical_time, label = 'carbon removed')
plt.plot(basal_oxygen, physical_time, label = 'oxygen adsorbed on basal plane')
plt.plot(defect_edge_oxygen, physical_time, label = 'oxygen adsrobed on defect edge')
plt.title('Simulation at %dK & %dPa'%(T_g,P_stag))
plt.xlabel('Physical time [ms]')
plt.ylabel('Number of atoms') 
plt.savefig("./Results/Plot_%dK_%dPa.png"%(T_g, P_stag))

