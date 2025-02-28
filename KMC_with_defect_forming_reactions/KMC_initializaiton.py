import numpy as np
import random
import math
import time
import itertools
from itertools import permutations  
import pandas as pd
import matplotlib.pyplot as plt

start_time = time.time()

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def within_4_neighs(site):
	all_C = []
	all_C_1 = []
	all_C_2 = []

	list_=[]
	for c in near_neigh[site]:
		all_C.extend(near_neigh[c])

	for c in all_C:
		all_C_1.extend(near_neigh[c])
	for c in all_C_1:
		all_C_2.extend(near_neigh[c])

	[list_.append(x) for x in all_C_2+all_C_1+all_C+near_neigh[site] if x not in list_]
	return list_

def periodic(site_a,site_b,coord_list,max_x,max_y):
	dist_x = abs(coord_list[site_a][0]-coord_list[site_b][0])
	dist_y = abs(coord_list[site_a][1]-coord_list[site_b][1])
	#lay = int(coord_list[site_a][3])-1
	#max_x = max_x_list[lay]
	#max_y = max_y_list[lay]
	#print(dist_x)
	#print(dist_y)
	if dist_x <(max_x/2) and dist_y <(max_y/2):
		return np.sqrt((coord_list[site_a][0]-coord_list[site_b][0])**2 + (coord_list[site_a][1]-coord_list[site_b][1])**2)
	elif dist_x <(max_x/2) and dist_y >(max_y/2):
		return np.sqrt((coord_list[site_a][0]-coord_list[site_b][0])**2 + (abs(coord_list[site_a][1]-coord_list[site_b][1])-max_y)**2)
	elif dist_x >(max_x/2) and dist_y <(max_y/2):
		return np.sqrt((abs(coord_list[site_a][0]-coord_list[site_b][0])-max_x)**2 + (coord_list[site_a][1]-coord_list[site_b][1])**2)
	else:
		return np.sqrt((abs(coord_list[site_a][0]-coord_list[site_b][0])-max_x)**2 + (abs(coord_list[site_a][1]-coord_list[site_b][1])-max_y)**2)

def isclose(a, b, tol):
	if abs(a-b) <= tol:
	    return True
	else:
		return False

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
				#print("Other Option")
		elif len(list_n) == 1:
			return ["SBC",list_n] #SBC - Single bonded carbon
		else:
			return ["Island",list_n]
	else:
		return ["Not in edge",list_n]

def armchair_pair(site):
	return next(p for p in near_neigh[site] if 'dg' not in flag[p] and [flag[k] for k in near_neigh[p]].count(['dg']) == 1)


def add_upper_apron_x(Upper_apron_x, Upper_apron_y, Lower_apron_x, Lower_apron_y, coordinates_per):
	for l in layers:
		for i in Upper_apron_x[l]:
			"""
			if (l % 2) == 0:
				coordinates_per=coordinates_per+[[coordinates[i][0]-(Upper_limit_x[l]+0.5+1.42),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
			else:
			"""
			coordinates_per=coordinates_per+[[coordinates[i][0]-(max_C_x+0.71),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Upper_apron_y[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0],coordinates[i][1]-(max_C_y+1.22976),coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Lower_apron_x[l]:
			"""
			if (l % 2) == 0:
				coordinates_per=coordinates_per+[[coordinates[i][0]+(Upper_limit_x[l]+0.5+1.42),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
			else:
			"""
			coordinates_per=coordinates_per+[[coordinates[i][0]+(max_C_x+0.71),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Lower_apron_y[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0],coordinates[i][1]+(max_C_y+1.22976),coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
	return coordinates_per

# parameters

size = 10
no_of_layers = 1

print("For size %d"%size)
print("For number of layers %d"%no_of_layers)
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
coordinates = np.copy(np.hstack((coordC_x, coordC_y, coordC_z, layernno, id_no)))
coordinates = np.stack( sorted(coordinates, key=lambda x: -x[2]), axis=0 )
for k in range(len(coordinates)):
	coordinates[k][4] = k
sites = len(coordinates)
no_C_each_layer = int(sites/no_of_layers)
layers = range(no_of_layers)


coordinates_O_2 = np.copy(np.hstack((coordC_x, coordC_y, coordC_z+1.275, layernno)))
coordinates_O_2 = np.stack( sorted(coordinates_O_2, key=lambda x: -x[2]), axis=0 )

max_C_x = max(coordC_x)[0]
max_C_y = max(coordC_y)[0]

#Oxygen lattice
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

coord_x=np.array(np.reshape(coord_x, (len(coord_x),1)))
coord_y=np.array(np.reshape(coord_y, (len(coord_y),1)))
coord_z=np.array(np.reshape(coord_z, (len(coord_z),1)))
layernno_O = np.array([no_of_layers-1-int(math.ceil((x-1.275)/3.4)) for x in coord_z]) 
layernno_O = np.array(np.reshape(layernno_O, (len(layernno_O),1)))

coordinates_O_1 = np.copy(np.hstack((coord_x, coord_y, coord_z, layernno_O)))
coordinates_O_1 = np.stack(sorted(coordinates_O_1, key=lambda x: -x[2]), axis=0 )
coordinates_O_3 = np.copy(coordinates)

sites_O_1 = len(coordinates_O_1) #Epoxy
sites_O_2 = len(coordinates_O_2) #Adatom
sites_O_3 = len(coordinates_O_3) #In-plane
no_ep_each_layer = int(sites_O_1/no_of_layers)
# Max Lengths

max_O_x = max(coord_x)[0]
max_O_y = max(coord_y)[0]

size_x =  (3*1.42*size*2)
size_y = (np.sqrt(3))*1.42*size*3
size_z = no_of_layers*3.4

full_size_c = size_x*size_y*size_z
unit_size =  3*1.42*(np.sqrt(3))*1.42 * 3.4

size_x =  (3*1.42*size*2)
size_y = (np.sqrt(3))*1.42*size*3
size_z = no_of_layers*3.4

no_of_catoms = int(full_size_c*4/unit_size)
no_of_oatoms = (no_of_catoms/4)*6

"""
f1=open("bonds.reaxc","r")
file=f1.readlines()

near_neigh=[0]*no_of_catoms
epoxy_O=[0]*no_of_oatoms
epoxy_O_move=[0]*no_of_oatoms

count=0
for line in file:
    count += 1
linesno=count

for line in file[7:count]:
    idno = int(line.split()[0])-1
    type_= int(line.split()[1])
    bonds=int(line.split()[2])
    if type_ == 1:
        for k in range(bonds):
            a = int(line.split()[3+k])
            if a<no_of_catoms:
                near_neigh[idno]=a
    else:
        idno = idno-no_of_catoms
        for k in range(bonds):
            a = int(line.split()[3+k])
            if a<no_of_catoms:
                epoxy_O[idno]=a
            else:
                epoxy_O_move[idno]=a
        
f1.close()
"""
Upper_limit_x = [0]*no_of_layers
Upper_limit_y = [0]*no_of_layers
Lower_limit_x = [0]*no_of_layers
Lower_limit_y = [0]*no_of_layers

Upper_apron_x = [0]*no_of_layers
Upper_apron_y = [0]*no_of_layers
Lower_apron_x = [0]*no_of_layers
Lower_apron_y = [0]*no_of_layers

for l in layers:
	Upper_limit_x[l] = max([coordinates[x][0] for x in range(int((l/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites)))])
	Upper_limit_y[l] = max([coordinates[x][1] for x in range(int((l/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites)))])
	Lower_limit_x[l] = min([coordinates[x][0] for x in range(int((l/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites)))])
	Lower_limit_y[l] = min([coordinates[x][1] for x in range(int((l/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites)))])
	#upper_apron_x =[ coordinates[x][0]>Upper_limit_x for x in range(int(((l)/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites)))]
	Upper_apron_x[l] = [ x for x in range(int(((l)/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites))) if coordinates[x][0]>Upper_limit_x[l]-0.5 and coordinates[x][3] == 1.0*l]
	Upper_apron_y[l] = [ x for x in range(int(((l)/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites))) if coordinates[x][1]>Upper_limit_y[l]-0.5 and coordinates[x][3] == 1.0*l]
	Lower_apron_x[l] = [ x for x in range(int(((l)/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites))) if coordinates[x][0]<Lower_limit_x[l]+0.5 and coordinates[x][3] == 1.0*l]
	Lower_apron_y[l] = [ x for x in range(int(((l)/no_of_layers)*(sites)),int(((l+1)/no_of_layers)*(sites))) if coordinates[x][1]<Lower_limit_y[l]+0.5 and coordinates[x][3] == 1.0*l]


# LIST 1
coordinates_per = np.copy(coordinates)
coordinates_per= coordinates_per.tolist()
coordinates_per = add_upper_apron_x(Upper_apron_x, Upper_apron_y, Lower_apron_x, Lower_apron_y, coordinates_per)
coordinates_per= np.array(coordinates_per)
coordinates_per = np.stack(sorted(coordinates_per, key=lambda x: -x[2]), axis=0 )
sites_per = len(coordinates_per)
sites_per_ind = [int(coordinates_per[x][4]) for x in range(0,sites_per)]

# LIST 1
# Near Neighbor list
near_neigh=[]
sites = len(coordinates)
for i in range(sites):
	listap = [p for p in range(sites) if (coordinates[p][1]-coordinates[i][1])**2 + (coordinates[p][2]-coordinates[i][2])**2 + coordinates[p][3]==coordinates[i][3] <= 1.43**2 and p!=i]
	near_neigh.append(listap)

# Lower Layer Exposure

print("Done")

check_above = [[]]*no_C_each_layer
for l in layers[1:no_of_layers]:
	#print("%d | Near Neighbors"%l)
	sites_above = (list(range(sites_per)))[0:int((l/no_of_layers)*(sites_per))]
	#sites_above = (list(range(sites_per)))[int(((l-1)/no_of_layers)*(sites_per)):int((l/no_of_layers)*(sites_per))]
	sites_check = (list(range(sites)))[int((l/no_of_layers)*(sites)):int(((l+1)/no_of_layers)*(sites))]
	for x in sites_check:
		check_above.append([sites_per_ind[w] for w in sites_above if (coordinates_per[w][0]-coordinates[x][0])**2 + (coordinates_per[w][1]-coordinates[x][1])**2 <= 1.43**2])
		#print("%d | Check Above"%x)

print("Done")

# LIST 3
check_below=[]
for l in layers[1:no_of_layers]:
	sites_check = (list(range(sites)))[int(((l-1)/no_of_layers)*(sites)):int((l/no_of_layers)*(sites))]
	sites_below = (list(range(sites_per)))[int((l/no_of_layers)*(sites_per)):(sites_per)]
	for x in sites_check:
		check_below.append([sites_per_ind[w] for w in sites_below if (coordinates_per[w][0]-coordinates[x][0])**2 + (coordinates_per[w][1]-coordinates[x][1])**2 <= 1.43**2])
		#print("%d | Check Below"%x)

check_below.extend([[]]*no_C_each_layer)

print("Done")

# LIST 4 
#C-C-O for epoxy
epoxy_O=[]
for i in range(sites_O_1):
	#print("%d | Epoxy O |"%i)
	listap = [sites_per_ind[p] for p in range(sites_per) if ((coordinates_per[p][0]-coordinates_O_1[i][0])**2 + (coordinates_per[p][1]-coordinates_O_1[i][1])**2 <= 0.72**2 and coordinates_per[p][3] == coordinates_O_1[i][3])]
	listap.sort()
	epoxy_O.append(listap)

edge_epoxy =[]
for i in range(sites_O_1):
	if len(epoxy_O[i])<2:
		edge_epoxy.append(i)

epoxy_O = [x for x in epoxy_O if epoxy_O.index(x) not in edge_epoxy]
coordinates_O_1= np.delete(coordinates_O_1,edge_epoxy,axis=0)
sites_O_1 = len(coordinates_O_1) #Epoxy

epoxy_O=[]
for i in range(sites_O_1):
    listap = [sites_per_ind[p] for p in range(sites_per) if ((coordinates_per[p][0]-coordinates_O_1[i][0])**2 + (coordinates_per[p][1]-coordinates_O_1[i][1])**2 <= 0.72**2 and coordinates_per[p][3] == coordinates_O_1[i][3])]
    listap.sort()
    epoxy_O.append(listap)

print("Done")

Arm_chair_neigh =[]
for site in range(sites):
	#print("%d | Arm Chair |"%site)
	all_C = []
	all_C_1 = []

	for c in near_neigh[site]:
		all_C.extend(near_neigh[c])
	all_C =list(set(all_C))
	for c in all_C:
		all_C_1.extend(near_neigh[c])
	Arm_chair_neigh.append(list(set([p for p in all_C_1 if all_C_1.count(p)==2])))

print("Done")

epoxy_list = []
ZigZag_neigh = []
list_=[]

near_neigh_C_O_7 =[]
near_neigh_C_8 =[]
for site in range(sites):
	#print("%d | ZigZag | Near Neigh C O 5 | Near Neigh C 5"%site)
	all_C = []
	all_C_1 = []
	all_C_2 = []
	all_C_3 = []
	all_C_4 = []
	all_C_5 = []
	all_C_6 = []
	list_=[]
	all_ep =[]

	for c in near_neigh[site]:
		all_C.extend(near_neigh[c])
	all_C =list(set(all_C))
	for c in all_C:
		all_C_1.extend(near_neigh[c])
	all_C_1 =list(set(all_C_1))
	for c in all_C_1:
		all_C_2.extend(near_neigh[c])
	all_C_2 =list(set(all_C_2))
	for c in all_C_2:
		all_C_3.extend(near_neigh[c])
	all_C_3 =list(set(all_C_3))
	for c in all_C_3:
		all_C_4.extend(near_neigh[c])
	all_C_4 =list(set(all_C_4))
	for c in all_C_4:
		all_C_5.extend(near_neigh[c])
	all_C_5 =list(set(all_C_5))
	for c in all_C_5:
		all_C_6.extend(near_neigh[c])
	all_C_6 =list(set(all_C_6))
	
	[list_.append(x) for x in all_C_4+all_C_3+all_C_2+all_C_1+all_C+near_neigh[site]+[site] if x not in list_]
	[[all_ep.append(x) for x in [epoxy_list[k][0][1],epoxy_list[k][1][1],epoxy_list[k][2][1]] if x not in all_ep] for k in list_]
	near_neigh_C_O_7.append(all_ep)
	[list_.append(x) for x in all_C_6+all_C_5 if x not in list_]
	near_neigh_C_8.append(list_)

	ZigZag_neigh.append([x for x in all_C if x != site])
	list_.append(list(set(all_C_3+all_C_2+all_C_1+all_C+near_neigh[site]+[site])))
	epoxy_list.append([[near_neigh[site][0],epoxy_O.index(sorted((near_neigh[site][0],site)))],[near_neigh[site][1],epoxy_O.index(sorted((near_neigh[site][1],site)))],[near_neigh[site][2],epoxy_O.index(sorted((near_neigh[site][2],site)))]])

print("Done")

lst_1 = list(zip(near_neigh, epoxy_list, check_above, check_below, ZigZag_neigh, Arm_chair_neigh, near_neigh_C_O_5, near_neigh_C_5))
df_1 = pd.DataFrame(lst_1, columns =['Near Neighbors','Epoxy List', 'Check Above', 'Check Below', 'ZigZiag Neighbors', 'ArmChair Neighbors', 'Near Neighbors C O 5','Near Neighbors C 5'], dtype = int)
df_1.to_pickle("Carbon_%dx%d_%d.pkl"%(size,size,no_of_layers))

print("Done")

# LIST 4 
#C-C-O for epoxy

epoxy_C_move = np.empty(sites_O_1, dtype=object) 
epoxy_O_move = np.empty(sites_O_1, dtype=object) 
epoxy_C_orient = np.empty(sites_O_1, dtype=object) 
epoxy_O_orient = np.empty(sites_O_1, dtype=object) 
near_neigh_8 =[]
near_neigh_O_7=[]

for site_O in range(sites_O_1):
	#print("%d | near_neigh_O_7 |"%site_O)
	all_C_ = [[],[]]
	all_C_1_ = [[],[]]
	all_C_2_ = [[],[]]
	all_C_3_ = [[],[]]
	all_C_4_ = [[],[]]
	all_C_5_ = [[],[]]
	all_C_6_ = [[],[]]

	all_C = [[],[]]
	all_C_1 = [[],[]]
	all_C_2 = [[],[]]
	all_C_3 = [[],[]]
	all_C_4 = [[],[]]
	all_C_5 = [[],[]]
	all_C_6 = [[],[]]

	list_=[]
	all_ep =[]
	
	for i in range(2):
		for c in near_neigh[epoxy_O[site_O][i]]:
			all_C_[i].extend(near_neigh[c])
		[all_C[i].append(x) for x in all_C_[i] if x not in all_C[i]]
		for c in all_C[i]:
			all_C_1_[i].extend(near_neigh[c])
		[all_C_1[i].append(x) for x in all_C_1_[i] if x not in all_C_1[i]]
		for c in all_C_1[i]:
			all_C_2_[i].extend(near_neigh[c])
		[all_C_2[i].append(x) for x in all_C_2_[i] if x not in all_C_2[i]]
		for c in all_C_2[i]:
			all_C_3_[i].extend(near_neigh[c])
		[all_C_3[i].append(x) for x in all_C_3_[i] if x not in all_C_3[i]]
		for c in all_C_3[i]:
			all_C_4_[i].extend(near_neigh[c])
		[all_C_4[i].append(x) for x in all_C_4_[i] if x not in all_C_4[i]]
		for c in all_C_4[i]:
			all_C_5_[i].extend(near_neigh[c])
		[all_C_5[i].append(x) for x in all_C_5_[i] if x not in all_C_5[i]]
		for c in all_C_5[i]:
			all_C_6_[i].extend(near_neigh[c])
		[all_C_6[i].append(x) for x in all_C_6_[i] if x not in all_C_6[i]]

	[list_.append(x) for x in all_C_5[1]+all_C_4[1]+all_C_3[1]+all_C_2[1]+all_C_1[1]+all_C[1]+all_C_5[0]+all_C_4[0]+all_C_3[0]+all_C_2[0]+all_C_1[0]+all_C[0]+near_neigh[epoxy_O[site_O][0]]+near_neigh[epoxy_O[site_O][1]]+epoxy_O[site_O] if x not in list_]
	[[all_ep.append(x) for x in [epoxy_list[k][0][1],epoxy_list[k][1][1],epoxy_list[k][2][1]] if x not in all_ep] for k in list_]
	near_neigh_O_7.append(all_ep)

	[list_.append(x) for x in all_C_6[1]+all_C_6[0] if x not in list_]
	for x in epoxy_O[site_O]:
		list_.remove(x)
	near_neigh_8.append(list_)

	epoxy_C_move[site_O] = [x for x in near_neigh[epoxy_O[site_O][0]]+near_neigh[epoxy_O[site_O][1]] if x not in epoxy_O[site_O]]
	epoxy_C_orient[site_O] =[x for x in all_C_1[0]+all_C_1[1] if x not in epoxy_O[site_O]]

	epoxy_O_move[site_O]=[epoxy_O.index(list(p)) for p in list(permutations(epoxy_C_move[site_O]+epoxy_O[site_O], 2)) if list(p) in epoxy_O and epoxy_O.index(list(p))!=site_O]
	epoxy_O_orient[site_O]=[epoxy_O.index(list(p)) for p in list(permutations(epoxy_C_orient[site_O], 2)) if list(p) in epoxy_O and epoxy_O.index(list(p))!=site_O]

print("Done")

lst_2 = list(zip(epoxy_O, epoxy_C_move, epoxy_O_move, epoxy_O_orient, near_neigh_6, near_neigh_O_6))

df_2 = pd.DataFrame(lst_2, columns =['Epoxy O', 'Epoxy C Move', 'Epoxy O Move', 'Epoxy O Orient', 'Near Neighbors 6','Near Neighbors O 6'], dtype = int)

df_2.to_pickle("Epoxies_%dx%d_%d.pkl"%(size,size,no_of_layers))

print("Done")
"""

output1 = pd.read_pickle("Carbon_%dx%d_%d.pkl"%(size,size,no_of_layers))
output2 = pd.read_pickle("Epoxies_%dx%d_%d.pkl"%(size,size,no_of_layers))
near_neigh = list(output1['Near Neighbors'])
epoxy_list = list(output1['Epoxy List'])
check_above = list(output1['Check Above'])
check_below = list(output1['Check Below'])
ZigZag_neigh = list(output1['ZigZiag Neighbors'])
Arm_chair_neigh = list(output1['ArmChair Neighbors'])
near_neigh_C_O_7 = list(output1['Near Neighbors C O 7'])
near_neigh_C_8 = list(output1['Near Neighbors C 8'])
epoxy_O = list(output2['Epoxy O'])
epoxy_C_move = list(output2['Epoxy C Move'])
epoxy_O_move = list(output2['Epoxy O Move'])
epoxy_O_orient = list(output2['Epoxy O Orient'])
epoxy_O_split = list(output2['Epoxy O Split'])
near_neigh_8 = list(output2['Near Neighbors 8'])
near_neigh_O_7 = list(output2['Near Neighbors O 7'])


"""

