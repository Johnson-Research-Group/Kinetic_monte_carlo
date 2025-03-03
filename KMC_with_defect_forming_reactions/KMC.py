import numpy as np
import random
import math
import time
import itertools
from itertools import chain
import pandas as pd
import cantera as ct
import sys

start_time = time.time()

class color:
	"""
	This function defines colors for text while printing
	"""

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

	"""
	This function finds neighboring sites to 'site' within four near neighbors.

	Parameters:
	site (int): carbon site index

	Returns:
	list: list of indices of near neighbors
	"""

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
	"""
	This function computes the distance between site_a & site_b given the periodicity of the graphene sheets.

	Parameters:
	site_a (int): carbon or epoxy site index
	site_b (int): carbon or epoxy site index
	coord_list (array): epoxy or carbon coordinate array
	max_x (float): maximum x value in coord_list
	max_y (float): maximum y value in coord_list

	Returns:
	float: The distance between site_a & site_b given periodicity.
	"""

	dist_x = abs(coord_list[site_a][0]-coord_list[site_b][0])
	dist_y = abs(coord_list[site_a][1]-coord_list[site_b][1])
	if dist_x <(max_x/2) and dist_y <(max_y/2):
		return np.sqrt((coord_list[site_a][0]-coord_list[site_b][0])**2 + (coord_list[site_a][1]-coord_list[site_b][1])**2)
	elif dist_x <(max_x/2) and dist_y >(max_y/2):
		return np.sqrt((coord_list[site_a][0]-coord_list[site_b][0])**2 + (abs(coord_list[site_a][1]-coord_list[site_b][1])-max_y)**2)
	elif dist_x >(max_x/2) and dist_y <(max_y/2):
		return np.sqrt((abs(coord_list[site_a][0]-coord_list[site_b][0])-max_x)**2 + (coord_list[site_a][1]-coord_list[site_b][1])**2)
	else:
		return np.sqrt((abs(coord_list[site_a][0]-coord_list[site_b][0])-max_x)**2 + (abs(coord_list[site_a][1]-coord_list[site_b][1])-max_y)**2)

def isclose(a, b, tol):

	"""
	This function determines if a is close to b within a tolerence value.

	Parameters:
	a (int or float): value
	b (int or float): value

	Returns:
	int: 1 indicates True and 0 indicates False
	"""

	if abs(a-b) <= tol:
		return 1
	else:
		return 0

def what_edge_carbon_is_it(site):

	"""
	This function determines what edge type 'site' is based on neighboring sites (zigzag, armchair, double zigzag, triple zigzag, sibgle-bonded carbon, island site, or not in edge).

	Parameters:
	site (int): carbon site index

	Returns:
	string: "ZZ" - zigzag, "AC" - armchair, "DZZ" - double zigzag, "TZZ" - triple zigzag, "SBC" - single-bonded carbon, "Island" - island site, "Not in edge" - basal plane site
	"""

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

def armchair_pair(site):

	"""
	This function determines the index of the site that is the armchair pair to 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	int: index of armchair pair of 'site'
	"""

	return next(p for p in near_neigh[site] if 'dg' not in flag[p] and [flag[k] for k in near_neigh[p]].count(['dg']) == 1)

def remove(site):

	"""
	This function flags 'site' as removed through function 'Epoxy_Ad_C_delete', and sets an covered sites in graphene layers below to be available for oxygen adsorption through function 'Lower_layer'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	for i in near_neigh[site]:
		Epoxy_Ad_C_delete(i)
		Lower_layer(i)

def epoxy_rit(site_O):

	"""
	This function updates the number sites available for reactions that deal with epoxy sites after list n1 is updated for site_O.

	Parameters:
	site_O (int): epoxy site index

	Returns:
	None
	"""

	before_ep_arr = (np.copy(changing_av_ep[site_O]))
	changing_av_ep[site_O]=np.array(n1)
	diff = (changing_av_ep[site_O]-before_ep_arr)
	sur_av_ep[:] = sur_av_ep[:] + diff

def adatom_rit(site):

	"""
	This function updates the number sites available for reactions that deal with carbon sites after list n2 is updated for site.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	before_ad_arr = np.copy(changing_av_ad[site])
	changing_av_ad[site]=np.array(n2)
	diff = changing_av_ad[site]-before_ad_arr
	sur_av_ad[:] = sur_av_ad[:] + diff

# Values needed to determine if epoxies are in a certain configuration in function ' epoxy_curr_orient'
pos = [2.45953,2.12997,3.25363,3.68925] 

def epoxy_curr_orient(site):

	"""
	This function classifies the orientaiton of one epoxy 'site' with neighboring epoxies based on 5 classifications depicted as Orient 1, Orient 2,Orient 3, Orient 4 in article.

	Parameters:
	site (int): epoxy site index

	Returns:
	list1, list2: list1 is a 2d list, where each list within the list contains the site indices of the epoxies in respective orientations (first list corresponds to Orient 1, second list corresponds to Orient 2..) with 'site'. list2 is a 1d list that contains either 1 or 0, where 1 indicates that an epoxy exists in the respective orientation (first value corresponds to Orient 1, second value corresponds to Orient 2..) with 'site'.
	"""

	ind_epo =[[int(p) for p in epoxy_O_orient[site] if isclose(periodic(p,site,coordinates_O_1,max_O_x+1.065,max_O_y+0.6149), j, 0.1) == 1 and coordinates_O_1[p][3] == coordinates_O_1[site][3] and flag_O_1[p] == 'ep' and [flag[epoxy_O[p][0]],flag[epoxy_O[p][1]]]==[['ep'],['ep']]] for j in pos]
	ind_epo[0] = [x for x in ind_epo[0] if len([j for j in epoxy_O_move[site] if j in epoxy_O_move[x]])>0]
	Orient = [1 if len(x)>0 else 0 for x in ind_epo+[[]]]

	if ind_epo == [[],[],[],[]]:
		ind_epo.append([0])
		Orient[4] = 1
	else:
		ind_epo.append([])

	return ind_epo, Orient

def epoxy_switch_orient(site, switch_site):
	
	"""
	This function finds the orientation of 'site' with 'switch_site' (Orient1, Orient2, Orient3, Orient4).

	Parameters:
	site (int): epoxy site index
	switch_site (int): epoxy site index

	Returns:
	list: a value of 1 or 0 indicating whether the epoxies are in the orientation corresponding to the index value of the element, respectively.
	"""

	val = [0]*5
	value=periodic(switch_site,site,coordinates_O_1,max_O_x+1.065,max_O_y+0.6149)
	if coordinates_O_1[switch_site][3] == coordinates_O_1[site][3]:
		val[1:4] = [isclose(value,pos[x], 0.1) for x in range(1,len(pos))]

		if isclose(value, 2.45953, 0.1) == 1 and len([x for x in epoxy_O_move[site] if x in epoxy_O_move[switch_site]])>0: val[0] = 1

	if val == [0]*5:
		val[4] = 1

	return val


def lactone_ether_group(site_C, flags):

	"""
	This function finds the site indices of relevant neighboring carbon atoms associated with the lactone-ether/ether-lactone-ether group.

	Parameters:
	site_C (int): carbon site index
	flags (string): 'la-eth' or 'eth-la-eth'

	Returns:
	int, int, int, int, int, int: site indices of carbon atoms associated with the lactone-ether or ether-lactone-ether groups.
	"""
	
	common_C = next(p for p in near_neigh[site_C]+[site_C] if 'ad' in flag[p] and flags in flag[p])

	stray_side_C = next(p for p in near_neigh[common_C] if 'stray' in flag[p])
	ether_C = next(p for p in near_neigh[common_C] if 'ether-side' in flag[p])
	lact_side_C = next(p for p in near_neigh[common_C] if 'lact-side' in flag[p])
	ether_C_side_C = next(p for p in near_neigh[ether_C] if p!=common_C and isclose(periodic(p,stray_side_C,coordinates,max_C_x+1.42,max_C_y+1.22975), 2.84, 0.1) == 1 and 'ether-side-side' in flag[p])
	cyc_eth_C = next(p for p in near_neigh[ether_C] if 'cyc-eth-C' in flag[p])
	return common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C, cyc_eth_C

def Epoxy_Ad_C_delete(C_atom):

	"""
	This function deletes the carbon site 'C_atom' and any oxygen atoms associated with it in the epoxy and carbonyl lattice. It also classifies the immediate near neighbors of the carbon atom as an edge site.

	Parameters:
	C_atom (int): carbon site index

	Returns:
	None
	"""

	lay = int(coordinates[C_atom][3])

	flag[C_atom] = ['dg']
	edge_C[C_atom] = 'dg'
	flag_O_2[C_atom] = 'dg'
	before_ad_arr = np.copy(changing_av_ad[C_atom])
	changing_av_ad[C_atom]=np.zeros((1,N2))
		
	sur_av_ad[:] = sur_av_ad[:] -before_ad_arr
	
	epoxy_del = [x[1] for x in epoxy_list[C_atom]]
	for q in epoxy_del:
		flag_O_1[q] = 'dg'
		before_ep_arr = (np.copy(changing_av_ep[q]))
		changing_av_ep[q]=np.zeros((1,N1))
		
		sur_av_ep[:] = sur_av_ep[:]-before_ep_arr

	for p in near_neigh[C_atom]:
		edge_C[p] = 'edg'

def Lower_layer(site_C):

	"""
	This function sets appropriate carbon atoms in layers directly below the carbon site 'site_C' as available if all relevant carbon atoms in the layers above have been removed.

	Parameters:
	site_C (int): carbon site index

	Returns:
	None
	"""

	make_avail_O =[]
	ll_sites_check = check_below[site_C]
	if len(ll_sites_check) !=0:
		for x in ll_sites_check:
			if cov_C[x] == 'cov':
				check_flag=[0]
				for w in check_above[x]:
					if flag[w] != ['dg']:
						break
					else:
						check_flag.append(1)

				if np.sum(check_flag) == len(check_above[x]): 
					cov_C[x] = 'uncov'
					for p in near_neigh[x]: cov_C[p] == 'uncov'
					[make_avail_O.append(epoxy_list[x][iter_][1]) for iter_ in range(3)]
	
	for p in make_avail_O:
		cov_ep[p] = 'uncov'
		if [flag[epoxy_O[p][0]],flag[epoxy_O[p][1]]]==[['pr'],['pr']] and flag_O_1[p] == 'pr':
			before_ep_arr = np.copy(changing_av_ep[p])
			changing_av_ep[p][0]=1
		
			if any( flag_O_1[y] == 'pr' and [flag[epoxy_O[y][0]],flag[epoxy_O[y][1]]]==[['pr'],['pr']] and cov_ep[y] == 'uncov' for y in epoxy_O_split[p]) == True:
				changing_av_ep[p][1]=1
			
			diff = (changing_av_ep[p]-before_ep_arr)
			sur_av_ep[:] = sur_av_ep[:]+diff

def epoxy_diffusion_ep_1(site):

	"""
	This function performs diffusion of an epoxy group to a carbon site that is part of a cyclic-ether group.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	cyc_eth_O = next(p for p in near_neigh[site] if flag_O_3[p] == 'cyc-eth')
	carb_cyc = [x for x in near_neigh[cyc_eth_O] if 'cyc-eth' in flag[x]]

	flag_O_3[cyc_eth_O] = 'na'
	cyc_ethers.remove(cyc_eth_O)
	for x in carb_cyc: flag[x] = ['carb_O']; flag_O_2[x] = 'carb_O'

	epoxy_remove_C(pivot_C)

def epoxy_diffusion_ep_2(site):

	"""
	This function performs diffusion of an epoxy group to a carbon site that is part of a lactone-ether group.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	common_C,lact_side_C,ether_C,stray_side_C,ether_C_side_C,cyc_eth_C  = lactone_ether_group(chosen, 'la-eth')
	list_ether_C = epoxy_list[ether_C]
	ether_O = [x[1] for x in list_ether_C if x[0]==common_C][0]
	ether_side_O=[x[1] for x in list_ether_C if x[0]==ether_C_side_C][0]

	flag[common_C],flag[lact_side_C],flag[ether_C],flag[stray_side_C],flag[ether_C_side_C], flag[cyc_eth_C],flag_O_2[common_C],flag_O_1[ether_O],flag_O_1[ether_side_O] = [['eth-la-eth','ad'], ['eth-la-eth','lact-side'], ['eth-la-eth','ether-side'], ['eth-la-eth','stray'], ['eth-la-eth','ether-side-side'],['eth-la-eth','cyc-eth-C'],'eth-la-eth','eth-la-eth','eth-la-eth']

	if flag[pivot_C]!= ['eth-la-eth','ether-side-side']:
			epoxy_remove_C(pivot_C)

def epoxy_diffusion_ep_3(site):

	"""
	This function performs diffusion of an epoxy group to a carbon site that is part of a carbonyl group.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	flag_O_2[pivot_C]='carb_O'
	flag[pivot_C].remove('ep')
	flag[pivot_C].append('carb_O')

def epoxy_diffusion_ep_fast(site):

	"""
	This function performs diffusion of an epoxy group when all carbon near neighbors around it are pristine sites.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch_O_ep = [p for p in epoxy_O_move[site] if chosen in epoxy_O[p]][0]
	flag_O_1[site] = 'pr'
	flag_O_1[switch_O_ep] = 'ep'
	for i in epoxy_O[site]: flag[i] = ['pr']
	for i in epoxy_O[switch_O_ep]: flag[i] = ['ep']

def split_O2_adsorption(site):

	"""
	This function performs the adsorption of two epoxies on the basal plane from molecular oxygen.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	chosen_O = np.random.choice([x for x in epoxy_O_split[site] if (flag_O_1[x]=='pr' and [flag[epoxy_O[x][0]],flag[epoxy_O[x][1]]]==[['pr'],['pr']])])

	for i in epoxy_O[site]: epoxy_make_ep_C(i)
	flag_O_1[site] = 'ep'
	flag_O_1[chosen_O] = 'ep'

	for i in epoxy_O[chosen_O]: epoxy_make_ep_C(i)

def adsorption_ep(site):

	"""
	This function performs the adsorption of an epoxy onto pristine carbon atoms.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	for i in epoxy_O[site]: epoxy_make_ep_C(i)
	flag_O_1[site] = 'ep'

def Recomb_O2_Desorp_ep(site):

	"""
	This function performs the recombination of two epoxies into gaseous molecular oxygen.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch = int(np.random.choice((epoxy_curr_orient(site)[0])[1]))

	for j in [site, switch]:
		flag_O_1[j] = 'pr' 
		for i in epoxy_O[j]: epoxy_remove_C(i)

def adsorption_ad(site):

	"""
	This function performs the adsorption of an oxygen atom onto a carbon site 'site'. Changes will be made to the flags in the carbon lattice and the carbonyl lattice. Since the carbon and carbonyl lattice have the same structure, the indices of the carbon and oxygen lattice are the same in both.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	flag_O_2[site] = 'ad'
	flag[site] = ['ad']

def epoxy_remove_C(site_C):

	"""
	This function removes an epoxy flag from a carbon site 'site_C', without removing any other flags the carbon site may have.

	Parameters:
	site_C (int): carbon site index

	Returns:
	None
	"""

	if flag[site_C] == ['ep']:
			flag[site_C] = ['pr']
	else:
		flag[site_C].remove('ep')

def epoxy_make_ep_C(site_C):

	"""
	This function flags 'site_C' as a carbon site with an epoxy adsorbed on it without affecting any flags it might have.

	Parameters:
	site_C (int): carbon site index

	Returns:
	None
	"""

	if flag[site_C] == ['pr']:
			flag[site_C] = ['ep']
	else:
		flag[site_C].append('ep')

def epoxy_move(site, switch_O):

	"""
	This function moves an epoxy from site with index 'site' to 'switch_O'.

	Parameters:
	site (int): epoxy site index
	switch_O (int): epoxy site index

	Returns:
	None
	"""

	flag_O_1[site] = 'pr'
	for i in epoxy_O[site]:
		epoxy_remove_C(i)
	flag_O_1[switch_O]='ep'
	for i in epoxy_O[switch_O]:
		epoxy_make_ep_C(i)

def ep_diff_1_2(site):

	"""
	This function moves an epoxy from site with index 'site' to a random epoxy index from the list 'switch_O_1' defined prior to calling this function in order to make the epoxy 'site' go from being in Orientation 1 to Orientation 2 with a nearby epoxy.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch_O = np.random.choice(switch_O_1)
	epoxy_move(site, switch_O)

def ep_diff_2_3(site):

	"""
	This function moves an epoxy from site with index 'site' to a random epoxy index from the list 'switch_O_2' defined prior to calling this function in order to make the epoxy 'site' go from being in Orientation 2 to Orientation 3 with a nearby epoxy.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch_O = np.random.choice(switch_O_2)
	epoxy_move(site, switch_O)

def ep_diff_1_4(site):

	"""
	This function moves an epoxy from site with index 'site' to a random epoxy index from the list 'switch_O_3' defined prior to calling this function in order to make the epoxy 'site' go from being in Orientation 1 to Orientation 4 with a nearby epoxy.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch_O = np.random.choice(switch_O_3)
	epoxy_move(site, switch_O)

def pos1_to_Epoxy_Ether(site):

	"""
	This function converts an epoxy 'site' and its neighboring epoxy group in Orientation 1 to an epoxy-ether group

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	flag_O_1[site] = 'pr'
	for i in epoxy_O[site]:
				flag[i] = ['pr']

	flag_O_1[switch_O] = 'ep-eth'
	flag_O_1[stay_put] = 'ep-eth'
	for i in epoxy_O[switch_O]: flag[i] = ['ep-eth','epoxy']
	for i in epoxy_O[stay_put]: flag[i] = ['ep-eth','eth'] if i != common_C else ['ep-eth','common']
	epox_C = next(p for p in epoxy_O[switch_O] if p != common_C)

	flag[opp_C] = ['ep-eth', 'opp-ad-C']
	flag[ether_C_side_C] = ['ep-eth', 'ether-C-side-C']
	flag[cyc_eth_C] = ['ep-eth','cyc-eth-C']

def epoxy_ether(site):
	
	"""
	This function finds the site indices of relevant neighboring carbon atoms associated with the epoxy-ether group, given an oxygen site 'site' in the epoxy-ether group.

	Parameters:
	site_C (int): carbon site index

	Returns:
	int, int, int, int, int, int, int, int, int: site indices of carbon and oxygen atoms associated with the epoxy-ether group.
	"""

	for i in epoxy_O[site]:
		if 'common' in flag[i]:
			common_C = i 
		else:
			pass
	list_common = epoxy_list[common_C]
	ether_C,ether_O = [x for x in list_common if flag_O_1[x[1]]=='ep-eth' and 'eth' in flag[x[0]]][0]
	epox_C,epox_O = [x for x in list_common if flag_O_1[x[1]]=='ep-eth' and 'epoxy' in flag[x[0]]][0]
	opp_ad_c = [p for p in near_neigh[common_C] if p!=epox_C and p!=ether_C][0]
	ether_C_side_C = [p for p in near_neigh[ether_C] if p!=common_C and isclose(periodic(p,opp_ad_c,coordinates,max_C_x+1.42,max_C_y+1.22975), 2.84, 0.1) == 1][0]
	cyc_eth_C = next(p for p in near_neigh[ether_C] if p not in [common_C,ether_C_side_C])
	return common_C, ether_C, epox_C, ether_O, epox_O, opp_ad_c, ether_C_side_C, cyc_eth_C

def Epoxy_Ether_to_pos1(site):

	"""
	This function converts an epoxy-ether group to two epoxies in Orientation 1, given oxygen site 'site' associated with the epoxy-ether.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	switch_O = [p for p in epoxy_O_move[site] if epoxy_switch_orient(p, ether_O)[0] == 1 and any(flag[g]==['pr'] for g in epoxy_O[p])==True][0]

	#common_C,ether_C,epox_C,ether_O,epox_O,opp_ad_c,ether_C_side_C = epoxy_ether(site)
	flag_O_1[epox_O] = 'pr'

	flag[opp_ad_c] = ['pr']
	flag[ether_C_side_C] = ['pr']
	flag[cyc_eth_C] = ['pr']

	flag_O_1[ether_O] = 'ep'
	for i in epoxy_O[ether_O]:
				flag[i] = ['ep']
	flag_O_1[switch_O] = 'ep'
	for i in epoxy_O[switch_O]:
				flag[i] = ['ep']

def Ep_Eth__Ep_Eth_Ep(site): #Epoxy Ether to Lactone Ether

	"""
	This function converts an epoxy-ether group to a lactone-ether group, given oxygen site 'site' associated with the epoxy-ether.

	Parameters:
	site (int): epoxy site index

	Returns:
	None
	"""

	common_C,ether_C,epox_C,ether_O,epox_O,opp_ad_c,ether_C_side_C,cyc_eth_C = epoxy_ether(site)

	other_C_ep = next(p for p in near_neigh[ether_C]+near_neigh[epox_C]+near_neigh[opp_ad_c]+near_neigh[ether_C_side_C]+near_neigh[cyc_eth_C] if flag[p]==['ep'])

	list_other_ep = epoxy_list[other_C_ep]
	other_C_ep_other = [x[0] for x in list_other_ep if flag[x[0]]==['ep'] and flag_O_1[x[1]]=='ep'][0]
	stray_ep = [x[1] for x in list_other_ep if x[0]==other_C_ep_other][0]
	
	flag[ether_C] = ['la-eth', 'eth','ether-side']
	j=np.random.choice([ether_C, epox_C])

	if j == ether_C:
		flag[j] = ['la-eth','eth','ether-side','double-eth']
		flag[epox_C] = ['la-eth','pr-la','lact-side']
	else:
		flag[epox_C] = ['la-eth','eth','lact-side']

	flag_O_1[stray_ep]='pr'
	epoxy_remove_C(other_C_ep)
	epoxy_remove_C(other_C_ep_other)
	flag_O_1[epox_O] = 'pr'

	flag[common_C] = ['la-eth','ad']
	flag_O_1[ether_O] ='la-eth'
	flag_O_2[common_C] ='la-eth'
	flag[opp_ad_c] = ['la-eth','eth','stray']
	flag[ether_C_side_C] = ['la-eth','ether-side-side']
	flag[cyc_eth_C] = ['la-eth','cyc-eth-C']

def la_eth_CO(site): #Lactone Ether to carbonyl and cyclic ether (CO)

	"""
	This function converts a lactone-ether group to a carbonyl and cyclic-ether group, given a carbon site 'site' associated with the lactone-ether.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C,cyc_eth_C = lactone_ether_group(site, 'la-eth')
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag_O_2[ether_C]='carb_O'
	flag[ether_C]=['carb_O']
	flag_O_3[site]='cyc-eth'
	cyc_ethers.append(site)
	flag[lact_side_C] =['cyc-eth']
	flag[stray_side_C] =['cyc-eth']
	flag[cyc_eth_C] = ['pr']
	flag[ether_C_side_C] = ['pr']

def eth_la_eth_CO(site): # Ether lactone ether to lactone at vacancy and cyclic ether (CO)

	"""
	This function converts a lactone-ether group to a lactone group at a vacancy, and a cyclic-ether, given a carbon site 'site' associated with the lactone-ether. This reaction produces a CO molecule.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C,cyc_eth_C = lactone_ether_group(site, 'eth-la-eth')
	list_ether_C = epoxy_list[ether_C]
	ether_O = [x[1] for x in list_ether_C if x[0]==site][0]
	other_ether = [x[1] for x in list_ether_C if x[0]==ether_C_side_C][0]

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag_O_2[ether_C]='la-vac'
	flag[ether_C] =['la-vac','ad']
	flag_O_1[other_ether] = 'la-vac'
	flag[ether_C_side_C] = ['la-vac','eth']
	flag[lact_side_C] = ['cyc-eth']
	flag[stray_side_C] = ['cyc-eth']
	flag[cyc_eth_C] = ['la-vac', 'cyc-eth-C']
	flag_O_3[site] = 'cyc-eth'
	cyc_ethers.append(site)

def la_eth_CO2(site): # lactone ether to carbonyl (CO2)

	"""
	This function converts a lactone-ether group to a carbonyl group, given a carbon site 'site' associated with the lactone-ether. This reaction produces a CO2 molecule.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C,cyc_eth_C = lactone_ether_group(site, 'la-eth')
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag_O_2[ether_C]='carb_O'
	flag[ether_C]=['carb_O']
	flag[stray_side_C] =['pr']
	flag[lact_side_C]=['pr']
	flag[ether_C_side_C]=['pr']
	flag[cyc_eth_C]=['pr']
	flag_O_1[ether_O] = 'dg'

def eth_la_eth_CO2(site): # Ether lactone ether to carbonyl and cyclic ether (CO2)

	"""
	This function converts a lactone-ether group to a carbonyl group and cyclic-ether group, given a carbon site 'site' associated with the lactone-ether. This reaction produces a CO2 molecule.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C,cyc_eth_C = lactone_ether_group(site, 'eth-la-eth')
	other_ether = [x[1] for x in epoxy_list[ether_C] if x[0]==ether_C_side_C][0]
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag_O_2[ether_C]='carb_O'
	flag[ether_C] =['carb_O']
	flag[lact_side_C] = ['cyc-eth']
	flag[stray_side_C] = ['cyc-eth']
	flag_O_3[site] = 'cyc-eth'
	cyc_ethers.append(site)
	flag[ether_C_side_C] =['pr']
	flag_O_1[other_ether] = 'pr'
	flag[cyc_eth_C]=['pr']

def la_vac_CO(site): #lactone at vacancy to cyclic ether (CO)

	"""
	This function converts a lactone group to a lactone group at a vacancy, and cyclic-ether group, given a carbon site 'site' associated with the lactone. This reaction produces a CO molecule.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	ether_C = next(p for p in near_neigh[site] if 'eth' in flag[p] and 'la-vac' in flag[p])
	other_C = next(p for p in near_neigh[site] if p!=ether_C and 'dg' not in flag[p])
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag_O_3[site] ='cyc-eth'
	cyc_ethers.append(site)
	flag[other_C] = ['cyc-eth']
	flag[ether_C] = ['cyc-eth']

def la_vac_epoxy_CO(site): #lactone at vacancy to 2 carbonyls (CO)

	"""
	This function converts a lactone group at a vacancy (in the presence of an epoxy) to 2 carbonyl groups, given a carbon site 'site' associated with the lactone. This reaction produces a CO molecule.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	ether_C = next(p for p in near_neigh[site] if 'eth' in flag[p] and 'la-vac' in flag[p])
	cyc_C = next(p for p in near_neigh[site] if 'cyc-eth-C' in flag[p] and 'la-vac' in flag[p])
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]
	imagin_O = [x[1] for x in epoxy_list[site] if x[0]==cyc_C][0]
	epoxy_C = [p for p in epoxy_C_move[ether_O]+epoxy_C_move[imagin_O] if flag[p]==['ep']][0]

	epoxy = [x[1] for x in epoxy_list[epoxy_C] if flag_O_1[x[1]]=='ep'][0]
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

	flag_O_1[epoxy]='pr'
	for i in epoxy_O[epoxy]: flag[i] = ['pr']

	flag_O_2[cyc_C] = 'carb_O'
	flag_O_2[ether_C] = 'carb_O'
	flag[cyc_C] = ['carb_O']
	flag[ether_C] = ['carb_O']

def la_vac_CO2(site): #lactone at vacancy to pristine (CO2)

	"""
	This function converts a lactone group at a vacancy to CO2 molecule, given a carbon site 'site' associated with the lactone.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	ether_C = next(p for p in near_neigh[site] if 'eth' in flag[p] and 'la-vac' in flag[p])
	other_C = next(p for p in near_neigh[site] if p!=ether_C and 'dg' not in flag[p])
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

	flag[ether_C] =['pr']

def la_vac_epoxy_CO2(site): #lactone at vacancy to cycliic ether (CO2)

	"""
	This function converts a lactone group at a vacancy (in the presence of an epoxy) to 2 cyclic-ethers, given a carbon site 'site' associated with the lactone. This reaction produces CO2.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	ether_C = next(p for p in near_neigh[site] if 'eth' in flag[p] and 'la-vac' in flag[p])
	cyc_C = next(p for p in near_neigh[site] if 'cyc-eth-C' in flag[p] and 'la-vac' in flag[p])
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]
	imagin_O = [x[1] for x in epoxy_list[site] if x[0]==cyc_C][0]
	epoxy_C = [p for p in epoxy_C_move[ether_O]+epoxy_C_move[imagin_O] if flag[p]==['ep']][0]

	epoxy = [x[1] for x in epoxy_list[epoxy_C] if flag_O_1[x[1]]=='ep'][0]
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

	flag_O_1[epoxy]='pr'
	for i in epoxy_O[epoxy]: flag[i] = ['pr']

	flag_O_3[site] ='cyc-eth'
	cyc_ethers.append(site)
	flag[cyc_C] = ['cyc-eth']
	flag[ether_C] = ['cyc-eth']

def la_eth_to_eth_la_eth(site):

	"""
	This function converts a lactone-ether group to an ether-lactone-ether group, given a carbon site 'site' associated with the lactone-ether.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	common_C,lact_side_C,ether_C,stray_side_C,ether_C_side_C,cyc_eth_C  = lactone_ether_group(site, 'la-eth')
	carbon_ep = next(p for p in list(itertools.chain(*[near_neigh[k] for k in lactone_ether_group(site, 'la-eth')])) if 'ep' in flag[p])
	list_eCsC = epoxy_list[carbon_ep]
	Other_ep_C = [x[0] for x in epoxy_list[carbon_ep] if (flag[x[0]]==['ep'] and flag_O_1[x[1]]=='ep')][0]
	epox_O = [x[1] for x in epoxy_list[Other_ep_C] if x[0]==carbon_ep][0]
	ether_O = [x[1] for x in epoxy_list[ether_C] if x[0]==common_C][0]
	ether_side_O= [x[1] for x in epoxy_list[ether_C_side_C] if x[0]==ether_C][0]

	flag_O_1[epox_O] = 'pr'
	for i in epoxy_O[epox_O]:
		epoxy_remove_C(i)

	flag[common_C] = ['eth-la-eth','ad']
	flag[lact_side_C] = ['eth-la-eth','lact-side']
	flag[ether_C]= ['eth-la-eth','ether-side']
	flag[stray_side_C] = ['eth-la-eth','stray']
	flag[ether_C_side_C] = ['eth-la-eth','ether-side-side']
	flag[cyc_eth_C] = ['eth-la-eth','cyc-eth-C']
	flag_O_2[common_C] = 'eth-la-eth'
	flag_O_1[ether_O] = 'eth-la-eth'
	flag_O_1[ether_side_O] = 'eth-la-eth'

def O2_edge_ad(site):

	"""
	This function performs adsorption of molecular oxygen onto a carbon edge site 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	flag[site] = ['carb_O2']
	flag_O_2[site] = 'carb_O2'

def O_edge_ad(site):

	"""
	This function performs adsorption of atomic oxygen onto the carbon edge site 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'

def O2_split_diff_ZZ1(site):

	"""
	This function performs dissociation and diffusion of O2 adsorbed onto a zigzag carbon edge site, 'site', onto a neighboring zigzag site.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	chosen_C_site = next(i for i in ZigZag_neigh[site] if flag[i]==['pr'] and what_edge_carbon_is_it(i)[0]=="ZZ")
	flag[chosen_C_site] = ['carb_O']
	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[chosen_C_site] = 'carb_O'

def O2_split_diff_AC1(site):

	"""
	This function performs dissociation and diffusion of O2 adsorbed onto a armchair carbon edge site, 'site', onto the pristine armchair pair site.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	other_arm_chair = armchair_pair(site)
	flag[site] = ['carb_O']
	flag[other_arm_chair] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[other_arm_chair] = 'carb_O'

def O2_split_diff_AC2(site):

	"""
	This function performs dissociation and diffusion of O2 adsorbed onto a armchair carbon edge site, 'site', onto a neighboring pristine armchair site that is not the armchair pair of 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	other_arm_chair = next(i for i in Arm_chair_neigh[site] if flag[i]==['pr'] and what_edge_carbon_is_it(i)[0]=="AC")
	flag[site] = ['carb_O']
	flag[other_arm_chair] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[other_arm_chair] = 'carb_O'

def O2_split_diff_ZZ2(site):

	"""
	This function performs dissociation and diffusion of O2 adsorbed onto a zigzag carbon edge site, 'site', as a carbonyl group, and an epoxy group.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	chosen_C = next(i for i in near_neigh[site] if flag[i]==['pr'] and any(flag[k]==['pr'] and edge_C[k] !='edg' for k in near_neigh[i])==True)
	ep_C2 = next(i for i in near_neigh[chosen_C] if flag[i]==['pr'] and edge_C[i] !='edg')

	ep_O = [x[1] for x in epoxy_list[ep_C2] if x[0]==chosen_C][0]
	flag_O_1[ep_O] ='ep'
	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	for i in epoxy_O[ep_O]: epoxy_make_ep_C(i)

def O2_des(site):

	"""
	This function performs the desoption of adsorbed molecular oxygen on defect edge site 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	flag[site] =['pr']
	flag_O_2[site] = 'pr'

def no_ep_ZZ_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a zigzag carbon edge site, 'site', in the presence of no epoxy groups.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def one_ep_ZZ_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a zigzag carbon edge site, 'site', in the presence of one epoxy group.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	chosen_C=next(p for p in near_neigh[site] if flag[p]==['ep'])
	list_chosen = epoxy_list[chosen_C]
	Other_ep_C = [k for k in near_neigh[chosen_C] if flag[k]==['ep']][0]
	epox_O = [x[1] for x in epoxy_list[Other_ep_C] if x[0]==chosen_C][0]

	other_epox_C= [x[0] for x in list_chosen if flag_O_1[x[1]]=='ep'][0]
	epox_O = [x[1] for x in list_chosen if x[0]==other_epox_C][0]
	flag_O_1[epox_O] = 'pr'
	for i in epoxy_O[epox_O]:
		epoxy_remove_C(i)

	flag[chosen_C] = ['carb_O']
	flag_O_2[chosen_C]='carb_O'

def two_ep_ZZ_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a zigzag carbon edge site, 'site', in the presence of two epoxy groups.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	for k in list_C:
		list_k= epoxy_list[k]
		other_epox_C= [x[0] for x in list_k if flag_O_1[x[1]]=='ep'][0]
		epox_O = [x[1] for x in list_k if x[0]==other_epox_C][0]

		flag_O_1[epox_O] = 'pr'
		for i in epoxy_O[epox_O]:
			epoxy_remove_C(i)

		flag[k] = ['carb_O']
		flag_O_2[k]='carb_O'

def one_AC_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a armchair carbon edge site, 'site', when the armchair pair site is pristine.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def two_AC_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a armchair carbon edge site, 'site', when the armchair pair site has a carbonyl group.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def one_AC_ep_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a armchair carbon edge site, 'site', when the armchair pair site is pristine, and there is one epoxy present in the vicinity of 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	epox_C = epox_C_1
	list_epox_C = epoxy_list[epox_C]
	other_epox_C= [x[0] for x in list_epox_C if flag_O_1[x[1]]=='ep'][0]
	epox_O = [x[1] for x in list_epox_C if x[0]==other_epox_C][0]
	flag_O_1[epox_O] = 'pr'

	for i in epoxy_O[epox_O]:
		epoxy_remove_C(i)

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	flag[epox_C] = ['carb_O']
	flag_O_2[epox_C]='carb_O'


def two_AC_ep_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a armchair carbon edge site, 'site', when the armchair pair site has a carbonyl group, and there is one epoxy present in the vicinity of 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	epox_C = next(p for p in [epox_C_1, epox_C_2] if flag[p] == ['ep'])
	list_epox_C = epoxy_list[epox_C]
	other_epox_C= [x[0] for x in list_epox_C if flag_O_1[x[1]]=='ep'][0]
	epox_O = [x[1] for x in list_epox_C if x[0]==other_epox_C][0]
	flag_O_1[epox_O] = 'pr'

	for i in epoxy_O[epox_O]:
		epoxy_remove_C(i)

	if site in near_neigh[epox_C]:
		Epoxy_Ad_C_delete(site)
		Lower_layer(site)
	else:
		arm_chair_pair = armchair_pair(site)
		Epoxy_Ad_C_delete(arm_chair_pair)
		Lower_layer(arm_chair_pair)
	flag[epox_C] = ['carb_O']
	flag_O_2[epox_C]='carb_O'

def one_DZZ(site):

	"""
	This function performs CO formation of a carbonyl group at a double zigzag site 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def two_DZZ(site):

	"""
	 This function performs CO formation of a carbonyl group at a double zigzag site 'site', where there are two carbonyl groups present at the double zigzag site.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def three_DZZ(site):

	"""
	This function performs CO formation of a carbonyl group at a double zigzag site 'site', where there are three carbonyl groups present at the double zigzag site.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def SBC_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a single-bonded carbon site 'site'.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def SBC_ep_CO(site):

	"""
	This function performs CO formation of a carbonyl group at a single-bonded carbon site 'site', in the presence of an epoxy group.

	Parameters:
	site (int): carbon site index

	Returns:
	None
	"""

	Epoxy_Ad_C_delete(site) # Epoxy group gets automatically deleted here
	Lower_layer(site)

# ============================================================================================
#                                        Parameters
# ============================================================================================

# Size unit of graphene sheets
size = int(sys.argv[1])

# No. of graphene layers
no_of_layers = int(sys.argv[2])

# KMC iterations to run
itr_final = int(sys.argv[3])

# KMC iterations to save
itr_save = int(sys.argv[4])

# Maximum walltime
walltime_max = int(sys.argv[5])

# Constants
k_b_ev = 8.6173 * 10**(-5) # Boltzmann's constant [eV K^-1]
k_b = 1.3806*(10**(-23)) # Boltzmann's constant [J/K]
h = 6.62607*(10**(-34)) # Planck's constant [Js]
m_O = 2.657*(10**(-26)) # Mass of oxygen atom [kg]
m_C = 1.99 * 10**(-26) # Mass of oxygen molecule [kg]
A_eff = 2.6199 * (10**-20) # Effective area for oxygen adsorption [m^2]

# Parameters
stick_coeff = 1 # Sticking coefficient of oxygen atoms for adsorption
T_s = T_g = 1500 # Temperature of surface and gas [K]
P_stag = 10000 # Pressure of gas mixture [Pa]

# Create cantera Solution object to represent a gas mixture gas1
gas1 = ct.Solution('gri30.yaml')

# Set temperature, pressure and species of gas mixture
gas1.TPX = T_g, P_stag, 'O2:0.21, N2:0.78, Ar:0.09'

# Equillibrate the gas mixture at set temperature and pressure
gas1.equilibrate('TP')

X = gas1['O','O2'].X # Mole fractions of atomic and molecular oxygen
P_O = X[0]*P_stag # Partial pressure of atomic oxygen
P_O2 = X[1]*P_stag # Partial pressure of molecular oxygen

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

# -------------------------------------- Print simulation details ---------------------------------

print("For size %d"%size)
print("For number of layers %d"%no_of_layers)
print("For temperature: %d"%T_s)
print("For Pressure: %d"%P_stag)
print("Partial Pressure O: %s"%str(P_O))
print("Partial Pressure O2: %s"%str(P_O2))

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

# Get number of reactions in 'Category 1' & 'Category 2'
N1 = len(Category1)
N2 = len(Category2)

# Create rate constant arrays for 'Category 1' & 'Category 2'
k1 = np.empty(N1, dtype=object) 
k2 = np.empty(N2, dtype=object) 

# ==================================================================================================
#                                      Reaction Rates [s^-1]
# ==================================================================================================

# -------------------------------------- Adsorption on basal plane --------------------------------
k1[0] = ((P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g))*stick_coeff # Adsorption of O on basal plane as epoxy
k1[1] = ((P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g))*np.exp(-1.7/(k_b_ev*T_s))*stick_coeff  # Dissociation of O2 & adsorption as two epoxies on the basal plane

# -------------------------------------- Desorption of epoxies ------------------------------------
k1[2] = (10**13)*np.exp(-1.130017/(k_b_ev*T_s)) # Recombination of two epoxies on the basal plane to desorb as O2

# -------------------------------------- Diffusion reactions --------------------------------------
k1[3]=k1[5] = (10**13)*np.exp(-0.737/(k_b_ev*T_s)) # Epoxy diffusion on basal plane
k1[4] =k1[6]= (10**13)*np.exp(-0.737/(k_b_ev*T_s)) # Epoxy diffusion between graphene sheets 
k1[7] = (10**13)*np.exp(-0.66/(k_b_ev*T_s)) # Epoxies shifting from Orientation 1 to Orientation 2
k1[8] = (10**13)*np.exp(-1.13/(k_b_ev*T_s)) # Epoxies shifting from Orientation 2 to Orientation 3
k1[9] = (10**13)*np.exp(-1.34/(k_b_ev*T_s)) # Epoxies shifting from Orientation 1 to Orientation 4
k1[10] = (10**13)*np.exp(-1.34/(k_b_ev*T_s)) # Epoxies shifting from Orientation 1 form an epoxy-ether group
k1[11] = (10**13)*np.exp(-0.95/(k_b_ev*T_s)) # Epoxy-ether group revert to two epoxies in Orientation 1
k1[12] = (10**13)*np.exp(-1.1/(k_b_ev*T_s)) # Lacton-ether group formation from epoxy-ether group

# -------------------------------------- Surface reactions ----------------------------------------
k2[0] = (k_b*T_s/h)*np.exp(-0.97/(k_b_ev*T_s)) # CO formation from lactone-ether group
k2[1] = (k_b*T_s/h)*np.exp(-0.5/(k_b_ev*T_s)) # CO formation from ether-latone-ether group
k2[2] = (k_b*T_s/h)*np.exp(-0.61/(k_b_ev*T_s)) # CO2 formation from lactone-ether group
k2[3] = (k_b*T_s/h)*np.exp(-0.6/(k_b_ev*T_s)) # CO2 formation from ether-lactone-ether group
k2[4] = (k_b*T_s/h)*np.exp(-2.98/(k_b_ev*T_s)) # CO formation from lactone at vacancy
k2[5] = (k_b*T_s/h)*np.exp(-1.93/(k_b_ev*T_s)) # CO formation from lactone at vacancy in the presence of a epoxy
k2[6] = (k_b*T_s/h)*np.exp(-2.46/(k_b_ev*T_s)) # CO2 formation from lactone at vacancy
k2[7] = (k_b*T_s/h)*np.exp(-0.58/(k_b_ev*T_s)) # CO2 from lactone at vacancy in the presence of a epoxy
k2[8] = (k_b*T_s/h) #Ether-lactone-ether formation from lactone-ether group in presence of epoxy

# -------------------------------------- Defect edge reactions ------------------------------------
k2[9] = (P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g)*stick_coeff  # O2 adsorption on zigzag edge site 
k2[10] = (P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g)*np.exp(-0.311/(k_b_ev*T_s))*stick_coeff # O2 adsorption on armchair edge site 
k2[11] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)*stick_coeff  # O adsorption on zigzag edge site 
k2[12] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)*stick_coeff # O adsorption on armchair edge site 
k2[13] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)*stick_coeff  # O adsorption on all other edge sites
k2[14] = 3.6*(10**12)*np.exp(-0.041/(k_b_ev*T_s)) # Dissociation & diffusion of O2 into two carbonyl groups adsorbed at neighboring zigzag sites
k2[15] = 1.5*(10**12)*np.exp(-0.694/(k_b_ev*T_s)) # Dissociation & diffusion of O2 into two carbonyl groups adsorbed at armchair pair sites
k2[16] = 3.5*(10**12)*np.exp(-0.238/(k_b_ev*T_s)) # Dissociation & diffusion of O2 into two carbonyl groups adsorbed at neighboring armchair sites
k2[17] = 1.2*(10**12)*np.exp(-1.513/(k_b_ev*T_s)) # Dissociation & diffusion of O2 into a carbonyl group adsorbed at a zigzag site and an epoxy on adjacent basal plane
k2[18] = 1.3*(10**14)*np.exp(-1.42/(k_b_ev*T_s)) # O2 desorption from zigzag site
k2[19] = 6.4*(10**14)*np.exp(-1.104/(k_b_ev*T_s)) # O2 desorption from armchair site
k2[20] = 1.2*(10**16)*np.exp(-3.622/(k_b_ev*T_s)) # CO formation from zigzag site, no epoxies present
k2[21] = 1.2*(10**16)*np.exp(-2.3/(k_b_ev*T_s)) # CO formation from zigzag site, one epoxy present
k2[22] = 1.2*(10**16)*np.exp(-2.3/(k_b_ev*T_s)) #  CO formation from zigzag site, two epoxies present
k2[23] = 1.0*(10**13)*np.exp(-2.743/(k_b_ev*T_s)) # CO formation from armchair site, pristine armchair pair
k2[24] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) # CO formation from armchair site, carbonyl present at armchair pair
k2[25] = 1.0*(10**13)*np.exp(-1.829/(k_b_ev*T_s)) # CO formation from armchair site, pristine armchair pair, one epoxy present
k2[26] = 1.0*(10**13)*np.exp(-1.136/(k_b_ev*T_s)) # CO formation from armchair site, carbonyl present at armchair pair, one epoxy present
k2[27] = 1.0*(10**13)*np.exp(-2.743/(k_b_ev*T_s)) # CO formation from carbonyl at double zigzag site
k2[28] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) # CO formation from carbonyl at double zigzag site, one additional carbonyl present
k2[29] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) # CO formation from carbonyl at double zigzag site, two additional carbonyl present
k2[30] = 1.0*(10**13)*np.exp(-2.17/(k_b_ev*T_s)) # CO formation from single-bonded carbon site
k2[31] = 1.0*(10**13)*np.exp(-0.447/(k_b_ev*T_s)) #CO formation from single-bonded carbon site, one epoxy present

# -------------------------------------- Island removal reactions ---------------------------------

# Set island removal rate for islands containing between 2 & 35 carbon atoms
for num in range(2,36):
	k2[num+31-1] = 1.0*(10**13)*np.exp(-(0.0751*num)/(k_b_ev*T_s)) 

# Set reation rate to floats
k1 = k1.astype('float')
k2 = k2.astype('float')
k_const = np.concatenate(((k1,k2)), axis=0)


# ==================================================================================================
#                                        Tracking sites
# ==================================================================================================

# ----------------------- Site status (surface group flag, exposure and edge) ----------------------

cov_C = ['uncov' if coordinates[i][3]==0 else 'cov' for i in range(sites)]
cov_ep = ['uncov' if coordinates_O_1[i][3]==0 else 'cov' for i in range(sites_O_1)]
edge_C = ['pr' for i in range(sites)]

flag = [['pr'] for i in range(sites)]
flag_O_1 = ['pr' for i in range(sites_O_1)]
flag_O_2 = ['pr' for i in range(sites_O_2)]
flag_O_3 = ['na' for i in range(sites_O_3)]

# ----------------------- Arrays to track available sites for each reaction ------------------------

changing_av_ad = np.zeros((sites_O_2,N2))
changing_av_ep = np.zeros((sites_O_1,N1))

sites_av_ad = [[] for i in range(N2)]
sites_av_ep = [[] for i in range(N1)]

# ==================================================================================================
#                                         Initialization
# ==================================================================================================

offset_t=0
offset_step=0
offset_full=0
site_O_before='na'
site_before='na'
switch_ad_before = [0]*no_of_layers
checker=0
cat=0
cyc_ethers = []

# Set all basal plane sites in top layer as available for epoxy adsorption
top_layer_O_1 = [i for i in range(sites) if coordinates_O_1[i][3]==0]

for i in top_layer_O_1:
	changing_av_ep[i][0]=1
	changing_av_ep[i][1]=1

sur_av_ep = np.array([len((np.where(changing_av_ep[:,n]==1)[0]).tolist()) for n in range(N1)])
sur_av_ad = np.array([len((np.where(changing_av_ad[:,n]==1)[0]).tolist()) for n in range(N2)])
sur_av = sur_av_ep.tolist()+sur_av_ad.tolist()

# Keep track of number of occurences of each reaction
counter_cat = [float(0)]*len(sur_av)


# ----------------------------- Re-run simulation from previous run --------------------------------

"""
# output3 = pd.read_pickle("Df_temp%d_19_7_22_%dPa_%dK_%s_%d_%d.pkl"%(set_,P_stag,T_g,str(size), no_of_layers, restart_timestep))
# flag, flag_O_1, flag_O_2, flag_O_3, offset_t, offset_step, offset_full, changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat = ((output3.tail(1)).values.tolist())[0]
"""

# ------------------------------ Initialize timestep and iteration ---------------------------------

t = offset_t # Physical time
time_step = int(offset_step) # KMC iteration


# ==================================================================================================
#                                       Begin KMC simualtion
# ==================================================================================================

# Run until no more reactions are available
while time_step<itr_final:

	time1 = time.time()

	# +++++++++++++++++++++++++++++++++++ [1] Find available sites +++++++++++++++++++++++++++++++++++

	# ----------------------------------- Reactions in Category 1 ------------------------------------

	# Identify epoxy sites within 7 neighbors of the previous epoxy site where the last reaction occurred.
	pool_ep=[]
	if site_O_before!='na': 
		pool_ep = near_neigh_O_7[site_O_before]+[site_O_before]
		before_ep_arr = (np.copy(changing_av_ep[site_O_before]))
		changing_av_ep[site_O_before]=np.zeros((1,N1))
		sur_av_ep = sur_av_ep-before_ep_arr

	if site_before!='na': 
		pool_ep= near_neigh_C_O_7[site_before]

	ind_epoxy = [x for x in pool_ep if flag_O_1[x]!='dg']
	ind_epoxy_pr = []
	ind_epoxy_ep = []
	ind_epoxy_eth_ep = [] 
	ind_epoxy_set = []

	# Organize sites into lists based on surface group status
	[ind_epoxy_pr.append(x) if (flag_O_1[x] == 'pr' and [flag[epoxy_O[x][0]],flag[epoxy_O[x][1]]] ==[['pr'],['pr']]) and cov_ep[x]=='uncov' else ind_epoxy_ep.append(x) if flag_O_1[x] == 'ep' else ind_epoxy_eth_ep.append(x) if (flag_O_1[x] == 'ep-eth' and any('epoxy' in flag[i] for i in epoxy_O[x])==True) else ind_epoxy_set.append(x) for x in ind_epoxy]

	n1 = [0]*N1
	[epoxy_rit(site_O) for site_O in ind_epoxy_set]

	ind_epoxy_pr_nor = []
	ind_epoxy_pr_split = []
	[ind_epoxy_pr_split.append(x) if any((flag_O_1[p] == 'pr' and [flag[epoxy_O[p][0]],flag[epoxy_O[p][1]]] ==[['pr'],['pr']]) and cov_ep[p]=='uncov' for p in epoxy_O_split[x])==True else ind_epoxy_pr_nor.append(x) for x in ind_epoxy_pr]

	n1[0]=1
	[epoxy_rit(site_O) for site_O in ind_epoxy_pr_nor]
	n1[1]=1
	[epoxy_rit(site_O) for site_O in ind_epoxy_pr_split]

	# Identify available reactions for epoxy sites within defined neighborhood
	for site_O in ind_epoxy_ep:

		n1 = [0]*N1  # List with values of 0 or 1, indicating whether a reaction in Category 1 is unavailable (0) or available (1), respectively.

		curr_list, curr_orient = epoxy_curr_orient(site_O)
		covered_val = cov_ep[site_O]

		n1[2]= 1 if curr_orient[1]==1 else 0
		C_move=[-1]*6

		if all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and all(flag[x]==['pr'] for x in near_neigh_8[site_O])==True: 
			# If all carbon sites are pristine around epoxy site

			n1[3] = 1 # and any(cov_C[x]=='uncov' for x in possible)==True # fast diffusion on exposed sheets
			#n1[4] = 1 if any(cov_C[x]=='cov' for x in epoxy_C_move[site_O])==True else 0 # fast diffusion between layers

		elif all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and all('ep' not in flag[x] for x in near_neigh_8[site_O])==True and all(flag[x]==['pr'] for x in near_neigh_8[site_O])!=True:
			# If there is atleast one non-pristine carbon site around the epoxy site

			occupied = [x for x in epoxy_C_move[site_O] if any(i in ['eth-la-eth','carb_O','carb_O2','la-vac','ep-eth','la-eth'] for i in flag[x])==True]
			n1[5] = 1 if len(occupied)<4 else 0 # and any(cov_C[x]=='uncov' for x in possible)==True # slow diffusion on exposed sheets
			#n1[6] = 1 if len(occupied)<4 and any(cov_C[x]=='cov' for x in possible)==True else 0 # slow diffusion between layers

		elif all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and any('ep' in flag[x] for x in near_neigh_8[site_O])==True:
			# If there are other epoxies present around the epoxy site

			# Check orientation of epoxy with other epoxies
			a,b,c,d=0,0,0,0
			if curr_orient[0]==1:
				for switch_O in epoxy_O_move[site_O]:
					pivot_C = next(p for p in epoxy_O[switch_O] if p not in epoxy_O[site_O])
					for p in curr_list[0]:
						list_orient = epoxy_switch_orient(int(p), switch_O)
						if list_orient[1] == 1 and flag[pivot_C]==['pr']:
							C_move[0] = epoxy_O[switch_O][0]
							C_move[1] = epoxy_O[switch_O][1]
							a=1
						if list_orient[3] == 1 and flag[pivot_C]==['pr']:
							C_move[2] = epoxy_O[switch_O][0]
							C_move[3] = epoxy_O[switch_O][1]
							c=1

				for stay_put in curr_list[0]:
					switch_O = [x for x in epoxy_O_move[site_O] if x in epoxy_O_move[stay_put]][0]
					common_C = [c for c in epoxy_O[switch_O] if c in epoxy_O[switch_O] and c in epoxy_O[stay_put]][0]
					opp_C = [p for p in near_neigh[common_C] if p not in epoxy_O[stay_put]+epoxy_O[switch_O]][0]
					ether_C = [p for p in epoxy_O[stay_put] if p != common_C][0]
					ether_C_side_C = next(p for p in near_neigh[ether_C] if p!=common_C and isclose(periodic(p,opp_C,coordinates,max_C_x+1.42,max_C_y+1.22975), 2.84, 0.1) == 1)
					cyc_eth_C = next(p for p in near_neigh[ether_C] if p not in [common_C,ether_C_side_C])
					if flag[opp_C]==['pr'] and flag[ether_C_side_C]==['pr'] and flag[cyc_eth_C]==['pr']:
						d=1

			if curr_orient[1]==1:
				for switch_O in epoxy_O_move[site_O]:
					pivot_C = next(p for p in epoxy_O[switch_O] if p not in epoxy_O[site_O])
					for p in curr_list[1]:
						list_orient = epoxy_switch_orient(int(p), switch_O)
						if epoxy_switch_orient(int(p), switch_O)[2] == 1 and flag[pivot_C]==['pr']:
								C_move[4] = epoxy_O[switch_O][0]
								C_move[5] = epoxy_O[switch_O][1]
								b=1

			C_move = [x for x in C_move if x not in epoxy_O[site_O] and x!=-1]
			occupied = [x for x in epoxy_C_move[site_O] if any(i in ['ep','eth-la-eth','carb_O','carb_O2','la-vac','ep-eth','la-eth'] for i in flag[x])==True]
			cyc_eth = []
			n1[5] = 1 if (len(C_move)+len(occupied))<4 else 0 # and any(cov_C[x]=='uncov' for x in possible)==True # slow diffusion on exposed sheets
			#n1[6] = 1 if (len(C_move)+len(occupied))<4 and any(cov_C[x]=='cov' for x in possible)==True else 0 # slow diffusion between layers

			n1[7] = 1 if a==1 else 0 # 1-Ep to 2-Ep 
			n1[8] = 1 if b==1 else 0 # 2-Ep to 3-Ep 
			n1[9] = 1 if c==1 else 0 # 1-Ep to 4-Ep 
			n1[10] = 1 if d==1 else 0 # 1-Ep to Epoxy-Ether

		else:
			pass

		epoxy_rit(site_O)

		# If the epoxy site is on island carbons with no neighbors, remove it 
		if all(flag[p] == ['dg'] for p in epoxy_C_move[site_O])==True:
			for i in epoxy_O[site_O]:
				Epoxy_Ad_C_delete(i)
				Lower_layer(i)

	# Identify available reactions for epoxy-ether sites within defined neighborhood
	for site_O in ind_epoxy_eth_ep:
		n1 = [0]*N1
		common_C,ether_C,epox_C,ether_O,epox_O,opp_ad_c,ether_C_side_C,cyc_eth_C = epoxy_ether(site_O)
		list_C_ = near_neigh[common_C]+near_neigh[ether_C]+near_neigh[epox_C]+near_neigh[opp_ad_c]+near_neigh[ether_C_side_C]+near_neigh[cyc_eth_C]
		list_C =[]
		[list_C.append(x) for x in list_C_ if x not in list_C]
		list_neigh = [x for x in list_C if x not in [common_C,ether_C,epox_C,opp_ad_c,ether_C_side_C,cyc_eth_C]]
		if all('dg' in flag[i] for i in list_neigh)==True:
			for p in [common_C,ether_C,epox_C,opp_ad_c,ether_C_side_C,cyc_eth_C]:
				Epoxy_Ad_C_delete(p)
				Lower_layer(p)
			
		else:
			common_C, ether_C, epox_C, ether_O, epox_O, opp_ad_c, ether_C_side_C,cyc_eth_C = epoxy_ether(site_O)
			n1[11] = 1 if len([1 for switch_O in epoxy_O_move[site_O] if any('ep-eth' in flag[x] for x in epoxy_O[site_O]) == True and epoxy_switch_orient(ether_O, switch_O)[0] == 1 and any(flag[p]==['pr'] for p in epoxy_O[switch_O])==True])>0 else 0 # Epoxy-Ether to 1-Ep
			if len([1 for k in near_neigh[ether_C]+near_neigh[epox_C]+near_neigh[opp_ad_c]+near_neigh[ether_C_side_C] +near_neigh[cyc_eth_C] if 'ep' in flag[k]])>0:
				n1[12] = 1 # Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation)

		epoxy_rit(site_O)

	# --------------------------------- Reactions in Category 2 --------------------------------------

	# Identify epoxy sites within 8 neighbors of the previous carbon site where the last reaction occurred.
	pool_ad=[]
	if site_O_before!='na': 
		pool_ad = near_neigh_8[site_O_before]+epoxy_O[site_O_before]
	if site_before!='na': 
		pool_ad = near_neigh_C_8[site_before]
		before_ad_arr = np.copy(changing_av_ad[site_before])
		changing_av_ad[site_before] = np.zeros((1,N2))
		
		sur_av_ad = sur_av_ad - before_ad_arr
	
	ind_ad_la = [x for x in pool_ad if flag_O_2[x] in ['la-vac','la-eth','eth-la-eth']]

	ind_ad_edg = [x for x in  pool_ad if edge_C[x] == 'edg' and flag[x]!=['dg'] and x not in ind_ad_la]
	edge_eps=[]
	[edge_eps.append(p) for p in ind_ad_edg if edge_C[p] =='edg' and flag[p] ==['ep']]
	ind_ad_edg = edge_eps+[p for p in ind_ad_edg if p not in edge_eps]
	n2 = [0]*N2

	# Identify available reactions for lactone sites within defined neighborhood
	for site in ind_ad_la:
		n2 = [0]*N2
		if 'la-eth' in flag[site] or 'eth-la-eth' in flag[site]:
			ether_C = next(p for p in near_neigh[site] if 'ether-side' in flag[p])			
		else:
			ether_C = next(p for p in near_neigh[site] if 'eth' and 'la-vac' in flag[p])
			cyc_C = next(p for p in near_neigh[site] if 'cyc-eth-C' and 'la-vac' in flag[p])
			ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]
			imagin_O = [x[1] for x in epoxy_list[site] if x[0]==cyc_C][0]
		
		if flag_O_2[site] == 'la-vac':
			if all(flag[p] == ['dg'] for p in epoxy_C_move[ether_O])!=True:
				lact_list = [p for p in near_neigh[site] if flag[p] not in [['dg'],['la-vac','eth']]]
				if [flag[p] for p in near_neigh[site]].count(['dg']) ==1 and len(lact_list)!=0:
					lact_C = lact_list[0]
					n2[4] = 1 if flag[lact_C]==['la-vac', 'cyc-eth-C'] else 0 #lactone CO formation at VACANCY
					n2[6] = 1 if flag[lact_C]==['la-vac', 'cyc-eth-C'] else 0 #lactone CO2 formation at VACANCY
					n2[5] = 1 if any(flag[p]==['ep'] for p in epoxy_C_move[ether_O]+epoxy_C_move[imagin_O])==True else 0 #lactone-epoxy CO formation at VACANCY
					n2[7] = 1 if any(flag[p]==['ep'] for p in epoxy_C_move[ether_O]+epoxy_C_move[imagin_O])==True else 0 #lactone-epoxy CO2 formation at VACANCY
				else:
					Epoxy_Ad_C_delete(site)
					Lower_layer(site)
					flag[ether_C] =['carb_O']
					flag_O_2[ether_C] ='carb_O'
				
			else:
				Epoxy_Ad_C_delete(site)
				Lower_layer(site)
				Epoxy_Ad_C_delete(ether_C)
				Lower_layer(ether_C)

		if flag_O_2[site] == 'la-eth':
			n2[0] = 1 #Lactone-Ether CO formatiom
			n2[2] = 1 #lactone-ether CO2 formation
			n2[8] = 1 if any('ep' in flag[p] for p in list(itertools.chain(*[near_neigh[k] for k in lactone_ether_group(site, 'la-eth')])))==True else 0 #Lactone-ether to ether-lactone-ether

		if flag_O_2[site] == 'eth-la-eth':
			n2[1] = 1 #Ether-latone-ether CO formation
			n2[3] = 1 #ether-lactone-ether CO2 formation

		adatom_rit(site)

	# Identify available reactions for carbonyl sites within defined neighborhood
	for site in ind_ad_edg:
		n2 = [0]*N2
		check_edge_C = 0
		
		edge_C_site_type,edge_C_list = what_edge_carbon_is_it(site) # get edge type

		# -------------------------------------- Islands ---------------------------------------------
		add_chain = [x for x in near_neigh[site] if flag[x]!=['dg']]
		chain = [site] + add_chain
		while len(add_chain)!=0:
			temp_add_chain =[]
			old_chain = chain
			for elem in add_chain:
				temp_add_chain += [x for x in near_neigh[elem] if flag[x]!=['dg'] and x not in temp_add_chain]
			add_chain = [x for x in temp_add_chain if flag[x]!=['cyc-eth'] and x not in chain]
			if len(chain+add_chain)>=37:
				chain = [0]*37
				break
			else:
				chain+=add_chain

		if len(chain)==1 and flag[site]!=['cyc-eth']: #edge_C_site_type == "Island":
			remove(site)
			check_edge_C=1
		elif len(chain) in range(2,36):
			n2[len(chain)+31-1]=1
		else:
			pass

		# ---------------------------------------------------------------------------------------------

		if any(i in ['la-vac','cyc-eth'] for i in flag[site])!=True:
			if flag[site] ==['pr'] and cov_C[site]=='uncov':
				if edge_C_site_type == "ZZ" or edge_C_site_type == "DZZ" or edge_C_site_type == "TZZ":
					n2[9] = 1 if all(flag[x]==['dg'] for x in check_above[site])==True else 0
					n2[11] = 1 if all(flag[x]==['dg'] for x in check_above[site])==True else 0

				elif edge_C_site_type == "AC":
					n2[10] = 1 if all(flag[x]==['dg'] for x in check_above[site])==True else 0
					n2[12] = 1 if all(flag[x]==['dg'] for x in check_above[site])==True else 0

				elif edge_C_site_type == "SBC":
					n2[13] = 1 if all(flag[x]==['dg'] for x in check_above[site])==True else 0

				else:
					pass

			elif flag[site] ==['carb_O2']:
				if edge_C_site_type == "ZZ":

					n2[14] = 1 if any(flag[i]==['pr'] and what_edge_carbon_is_it(i)[0] =="ZZ" for i in ZigZag_neigh[site])==True else 0
					value =[0]
					for i in near_neigh[site]:
						if flag[i]==['pr']:
							zz_nn = [k for k in near_neigh[i] if k!=site]
							d_val = 1 if any(flag[k]==['pr'] and edge_C[k] !='edg' for k in zz_nn)==True else 0
							value.append(d_val)

					n2[17] = 1 if sum(value)>0 else 0
					n2[18] = 1 

				elif edge_C_site_type == "DZZ":
					n2[18] = 1 

				elif edge_C_site_type == "TZZ":
					remove(site)

				elif edge_C_site_type == "AC":
					n2[15] = 1 if flag[armchair_pair(site)] == ['pr'] else 0
					n2[16] = 1 if any(flag[i]==['pr'] and what_edge_carbon_is_it(i)[0]=="AC" for i in Arm_chair_neigh[site])==True else 0
					n2[19] = 1

				elif edge_C_site_type == "SBC":
					flag[site] = ['carb_O']
					flag_O_2[site] = 'carb_O'

				else:
					pass 

			elif flag[site] ==['carb_O']:	
				if edge_C_site_type == "ZZ":

					list_C=[]
					for i in near_neigh[site]:
						if 'dg' not in flag[i]:
							list_C.append(i)

					n2[20] = 1 if [flag[p] for p in list_C].count(['ep'])==0 else 0
					n2[21] = 1 if [flag[p] for p in list_C].count(['ep'])==1 else 0 
					n2[22] = 1 if all(flag[p] == ['ep'] for p in list_C) == True else 0

				elif edge_C_site_type == "DZZ":
					list_C =[i for i in near_neigh[site] if 'dg' not in flag[i]]
					flag_C =[flag[i] for i in near_neigh[site] if 'dg' not in flag[i]]
					n2[27] = 1 if flag_C.count(['pr'])==2 else 0
					n2[28] = 1 if flag_C.count(['pr'])==1 else 0
					n2[29] = 1 if flag_C.count(['pr'])==0 else 0

				elif edge_C_site_type == "TZZ":
					remove(site)
					
				elif edge_C_site_type == "AC":
					ac_pair = armchair_pair(site)
					n2[23] = 1 if flag[ac_pair] == ['pr'] else 0
					n2[24] = 1 if flag[ac_pair] == ['carb_O'] else 0
					epox_C_1=next(i for i in near_neigh[site] if 'dg' not in flag[i] and i != ac_pair)
					epox_C_2=next(i for i in near_neigh[ac_pair] if 'dg' not in flag[i] and i != site)
					n2[25] = 1 if flag[epox_C_1] == ['ep'] and flag[ac_pair] == ['pr'] else 0
					n2[26] = 1 if any(flag[p] == ['ep'] for p in [epox_C_1, epox_C_2])==True and flag[ac_pair] == ['carb_O'] else 0

				elif edge_C_site_type == "SBC":
					link_C = next(i for i in near_neigh[site] if 'dg' not in flag[i])
					list_C =[i for i in near_neigh[link_C] if i != site]
					n2[30] = 1
					n2[31] = 1 if any(flag[i]==['ep'] for i in list_C)==True else 0
				else:
					pass
			else:
				pass

		if check_edge_C==1:
			n2 = [0]*N2

		adatom_rit(site)

	# List containing number of site available for each reaction, n
	sur_av = sur_av_ep.tolist()+sur_av_ad.tolist()

	# Tally reaction rates (k * n)
	rates = np.absolute(np.einsum('i,i->i',k_const,sur_av))
	R = np.einsum('i,i->',k_const,sur_av)
	
	# +++++++++++++++++++++++++++ [2] Choose reaction based on probabilities +++++++++++++++++++++++++

	if R != 0:

		# Probabilities determined from reaction rates
		prob = [i/R for i in rates]
		#print("Probability:")
		#print(["%.15f|%s"%(prob[i], Category[i]) if prob[i]!=0 else 0 for i in range(len(sur_av))])
		cat = int(np.random.choice(list(range(N1+N2)),p=prob))
		#print("Category is %s"%Category[cat])
		counter_cat[cat] += 1

	else:
		cat = -1
		#print("Category is Not Available")

	site_O_before= 'na'
	site_before= 'na'
	
	# --------------------- Choose epoxy site if chosen reaction in Category 1 -----------------------

	if cat<N1 and cat>1:
		site_O = np.random.choice(np.where(changing_av_ep[:,cat]==1)[0])
		site_O_before=site_O

		# Initialize values before perfoming reaction

		if flag_O_1[site_O] == 'ep':
			C_move=[-1]*6
			switch_O_1=[]
			switch_O_2=[]
			switch_O_3=[]

			
			curr_list, curr_orient = epoxy_curr_orient(site_O)
			if curr_orient[1]==1:
				for switch_O in epoxy_O_move[site_O]:
					pivot_C = next(p for p in epoxy_O[switch_O] if p not in epoxy_O[site_O])
					for p in curr_list[1]:
						if epoxy_switch_orient(int(p), switch_O)[2] == 1 and flag[pivot_C]==['pr']:
							C_move[4] = epoxy_O[switch_O][0]
							C_move[5] = epoxy_O[switch_O][1]
							switch_O_2.append(switch_O)

			if curr_orient[0]==1:
				for switch_O in epoxy_O_move[site_O]:
					pivot_C = next(p for p in epoxy_O[switch_O] if p not in epoxy_O[site_O])
					for p in curr_list[0]:
						list_orient = epoxy_switch_orient(int(p), switch_O)
						if list_orient[1] == 1 and flag[pivot_C]==['pr']:
							C_move[0] = epoxy_O[switch_O][0]
							C_move[1] = epoxy_O[switch_O][1]
							switch_O_1.append(switch_O)
						if list_orient[3] == 1 and flag[pivot_C]==['pr']:
							C_move[2] = epoxy_O[switch_O][0]
							C_move[3] = epoxy_O[switch_O][1]
							switch_O_3.append(switch_O)

				for stay_put_ in curr_list[0]:
					switch_O_ = [x for x in epoxy_O_move[site_O] if x in epoxy_O_move[stay_put_]][0]
					common_C_ = [c for c in epoxy_O[switch_O_] if c in epoxy_O[stay_put_]][0]
					opp_C_ = [p for p in near_neigh[common_C_] if p not in epoxy_O[stay_put_]+epoxy_O[switch_O_]][0]
					ether_C_ = [p for p in epoxy_O[stay_put_] if p != common_C_][0]
					ether_C_side_C_ = next(p for p in near_neigh[ether_C_] if p!=common_C_ and isclose(periodic(p,opp_C_,coordinates,max_C_x+1.42,max_C_y+1.22975), 2.84, 0.1) == 1)
					cyc_eth_C_ = next(p for p in near_neigh[ether_C_] if p not in [common_C_,ether_C_side_C_])

					if flag[opp_C_]==['pr'] and flag[ether_C_side_C_]==['pr'] and flag[cyc_eth_C_]==['pr']:
						stay_put = stay_put_
						switch_O = switch_O_
						common_C = common_C_
						opp_C = opp_C_
						ether_C = ether_C_
						ether_C_side_C = ether_C_side_C_
						cyc_eth_C = cyc_eth_C_

			C_move = [x for x in C_move if x not in epoxy_O[site_O] and x!=-1]

		elif flag_O_1[site_O] == 'ep-eth' and any('epoxy' in flag[i] for i in epoxy_O[site_O])==True :
			common_C,ether_C,epox_C,ether_O,epox_O,opp_ad_c,ether_C_side_C,cyc_eth_C = epoxy_ether(site_O)
		else:
			pass

	# ---------------------- Choose carbon site if chosen reaction in Category 2 ---------------------

	elif cat>(N1-1) and cat< len(sur_av)-1:
		site = np.random.choice(np.where(changing_av_ad[:,cat-N1]==1)[0])
		site_before=site
		
		# Initialize values before perfoming reaction

		if flag_O_2[site] == 'la-vac':
			if all(flag[p] == ['dg'] for p in epoxy_C_move[ether_O])!=True:
				lact_list = [p for p in near_neigh[site] if 'cyc-eth-C' and 'la-vac' in flag[p]]
				if [flag[p] for p in near_neigh[site]].count(['dg']) ==1 and len(lact_list)!=0:
					lact_C = lact_list[0]

		elif edge_C[site] =='edg' and any(i in ['la-vac','cyc-eth'] for i in flag[site])!=True:
			add_chain = [x for x in near_neigh[site] if flag[x]!=['dg']]
			chain = [site] + add_chain
			while len(add_chain)!=0 and len(chain)<=37 and len(add_chain)!=0:
				temp_add_chain =[]
				old_chain = chain
				for elem in add_chain:
					temp_add_chain += [x for x in near_neigh[elem] if flag[x]!=['dg'] and x not in temp_add_chain]
				add_chain = [x for x in temp_add_chain if flag[x]!=['cyc-eth'] and x not in chain]
				chain+=add_chain

			edge_C_type,edge_C_list = what_edge_carbon_is_it(site)

			if flag[site] ==['carb_O2']:
				if edge_C_type == "ZZ":
					zz_nn =[]
					list_C=[]
					for i in near_neigh[site]:
						if flag[i]==['pr']:
							list_C.append(i)
							zz_nn.append(near_neigh[i])
					zz_nn = list(itertools.chain(*zz_nn))
					zz_nn = [i for i in zz_nn if i != site]

				else:
					pass 

			elif flag[site] ==['carb_O']:
				if edge_C_type == "ZZ":
					list_C=[]
					for i in near_neigh[site]:
						if 'dg' not in flag[i]:
							list_C.append(i)

				elif edge_C_type == "DZZ":
					list_C =[i for i in near_neigh[site] if 'dg' not in flag[i]]
					flag_C =[flag[i] for i in near_neigh[site] if 'dg' not in flag[i]]
					
				elif edge_C_type == "AC":
					epox_C_1=next(i for i in near_neigh[site] if 'dg' not in flag[i] and i != armchair_pair(site))
					epox_C_2=next(i for i in near_neigh[armchair_pair(site)] if 'dg' not in flag[i] and i != site)

				elif edge_C_type == "SBC":
					link_C = next(i for i in near_neigh[site] if 'dg' not in flag[i])
					list_C =[i for i in near_neigh[link_C] if i != site]
				else:
					pass
			else:
				pass
		else:
			pass
	else:
		pass

	# +++++++++++++++++++++++++++++++++++ [3] Perform chosen reaction ++++++++++++++++++++++++++++++++

	if cat == 0: # Adsoprtion O
		site_O = np.random.choice(np.where(changing_av_ep[:,cat]==1)[0])
		site_O_before=site_O
		adsorption_ep(site_O)

	elif cat== 1: # O2 Split Adsorption
		site_O = np.random.choice(np.where(changing_av_ep[:,cat]==1)[0])
		site_O_before=site_O
		split_O2_adsorption(site_O)

	elif cat == 2: #O2 Recombination & Desorption
		Recomb_O2_Desorp_ep(site_O)

	elif cat == 3 or cat ==4: #Epoxy diffusion fast
		chosen = np.random.choice(epoxy_C_move[site_O])
		epoxy_diffusion_ep_fast(site_O)

	elif cat == 5 or cat==6: #Epoxy diffusion slow
		occupied = [x for x in epoxy_C_move[site_O] if any(i in flag[x] for i in ['ep','eth-la-eth','carb_O','carb_O2','la-vac','la-eth','ep-eth'])==True]
		options = [p for p in epoxy_C_move[site_O] if p not in C_move + occupied]
		chosen = np.random.choice(options)
		switch_O_ep = [p for p in epoxy_O_move[site_O] if chosen in epoxy_O[p]][0]
		pivot_C = next(p for p in epoxy_O[switch_O_ep] if p!= chosen)
		previous_C = next(p for p in epoxy_O[site_O] if p!= pivot_C)

		epoxy_remove_C(previous_C) 
		flag_O_1[site_O] = 'pr'

		if 'cyc-eth' in flag[chosen]:
			epoxy_diffusion_ep_1(site_O)
		elif flag[chosen] ==['dg']:
			epoxy_diffusion_ep_3(site_O)
		else:
			epoxy_make_ep_C(chosen)
			flag_O_1[switch_O_ep] = 'ep'

	elif cat == 7: # 1-ep to 2-ep
		ep_diff_1_2(site_O)

	elif cat == 8:# 2-ep to 3-ep
		ep_diff_2_3(site_O)

	elif cat == 9: # 1-ep to 4-ep
		ep_diff_1_4(site_O)

	elif cat == 10: # 1-ep to epoxy-ether
		pos1_to_Epoxy_Ether(site_O)

	elif cat == 11: # epoxy-ether to 1-ep
		Epoxy_Ether_to_pos1(site_O)

	elif cat == 12: # epoxy-ether to lactone-ether
		Ep_Eth__Ep_Eth_Ep(site_O)

	elif cat == 13: # lactone-ether CO formation
		la_eth_CO(site)

	elif cat == 14: # ether-latone-ether CO formation
		eth_la_eth_CO(site)

	elif cat ==15: # lactone-ether CO2 formation
		la_eth_CO2(site)

	elif cat == 16:# ether-lactone-ether CO2 formation
		eth_la_eth_CO2(site)

	elif cat == 17: # lactone CO formation at vacancy
		la_vac_CO(site)

	elif cat == 18: # lactone-epoxy CO formation at vacancy
		la_vac_epoxy_CO(site)

	elif cat == 19: # lactone CO2 formation at vacancy
		la_vac_CO2(site)

	elif cat == 20: # lactone-epoxy CO2 formation at vacancy
		la_vac_epoxy_CO2(site)

	elif cat == 21: # lactone-ether to ether-lactone-ether
		la_eth_to_eth_la_eth(site)

	if cat == 22 or cat == 23: # O2 zigzag & armchair adsorption 
		O2_edge_ad(site)

	elif cat ==24 or cat == 25 or cat == 26: # O zigzag adsorption & O armchair adsorption & O carbonyl adsorption
		O_edge_ad(site)

	elif cat == 27: # O2 split diffusion ZZ1
		O2_split_diff_ZZ1(site)

	elif cat == 28: # O2 split diffusion AC1
		O2_split_diff_AC1(site)

	elif cat == 29: # O2 split diffusion AC2
		O2_split_diff_AC2(site)

	elif cat == 30: # O2 split diffusion ZZ
		O2_split_diff_ZZ2(site)
		
	elif cat == 31: # O2 desorption ZZ
		O2_des(site)

	elif cat==32: # O2 desorption AC
		O2_des(site)

	elif cat==33: # No epoxy ZZ CO formation
		no_ep_ZZ_CO(site)

	elif cat==34: # 1 epoxy ZZ CO formation
		one_ep_ZZ_CO(site)

	elif cat==35: # 2 epoxy ZZ CO formation
		two_ep_ZZ_CO(site)

	elif cat==36: # 1 carbonyl AC CO formation
		one_AC_CO(site)

	elif cat==37: # 2 carbonyl AC CO formation
		two_AC_CO(site)

	elif cat==38: # 1 carbonyl epoxy AC CO formation
		one_AC_ep_CO(site)

	elif cat==39: # 2 carbonyl epoxy AC CO formation
		two_AC_ep_CO(site)

	elif cat==40: # 1 carbonyl DZZ CO formation
		one_DZZ(site)

	elif cat==41: # 2 carbonyl DZZ CO formation
		two_DZZ(site)

	elif cat==42: # 3 carbonyl DZZ CO formation
		three_DZZ(site)

	elif cat==43: # SBC CO formation
		SBC_CO(site)

	elif cat==44: # SBC epoxy CO formation
		SBC_ep_CO(site)

	elif cat>=45 and cat<(36+45): # Island removal
		for elem in chain:
			Epoxy_Ad_C_delete(elem)
			Lower_layer(elem)

	else: #ZigZag Scatter & Armchair Scatter & Carbonyl O Scatter
		pass

	# Remove cyclic-ethers if part of island
	for site in cyc_ethers:
		attach_C = [p for p in near_neigh[site] if 'cyc-eth' in flag[p]]

		if all([flag[p] for p in near_neigh[i]].count(['dg'])==2 for i in attach_C)==True:
			next_attach = []
			for k in attach_C:
				next_attach.append(next(p for p in near_neigh[k] if 'dg' not in flag[p]))

			if all([flag[p] for p in near_neigh[i]].count(['dg'])==2 for i in next_attach)==True:
				for k in attach_C:
					Epoxy_Ad_C_delete(k)
					Lower_layer(k)
				for i in next_attach:
					Epoxy_Ad_C_delete(i)
					Lower_layer(i)
				flag_O_3[site]='dg'
				cyc_ethers.remove(site)

		elif all([flag[p] for p in near_neigh[i]].count(['dg'])==3 for i in attach_C)==True:
			flag_O_3[site]='dg'
			cyc_ethers.remove(site)
			for k in attach_C: 
				Epoxy_Ad_C_delete(k)
				Lower_layer(k)
		else:
			pass

	# ------------------------------------- Save simulation data -------------------------------------

	if time_step%itr_save==0:
		# Save every 1000000 KMC iterations
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3, t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])
		df1_4.to_pickle("./Df_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))
		print(" ========= KMC Iteration = %d ==========" %time_step)
	
	if cat==-1:
		# Save at last iteration
		print("All carbon atoms are removed")
		print(" ========= Final KMC Iteration = %d ==========" %time_step)
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3, t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])
		df1_4.to_pickle("./Df_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))

	if (time.time() - start_time)/3600 >walltime_max and checker==0:
		# Save after a certain amount of simulation time has passed
		print("Maximum walltime limit exceeded")
		print(" ========= KMC Iteration = %d ==========" %time_step)
		checker=1
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3,t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])
		df1_4.to_pickle("./Df_restart_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))

	# +++++++++++++++++++++++++++++++++++++ [4] Update physical time +++++++++++++++++++++++++++++++++

	if R!=0:
		r1 = random.random()
		t += -(np.log(r1))/(R*1.0)
	else:
		pass
	
	# --------------------------- Print computational time per iteration -----------------------------

	#timestep_end = time.time()
	#print("Time for timestep: %s seconds"%str(timestep_end -time1))
	#print("Time Elapsed: %s mins"%str((time.time() - start_time)/60))
	#print(" ========= t_step = %d ==========" %time_step)

	time_step +=1


# ==================================================================================================
#                                       Print simulation data
# ==================================================================================================

print(color.BOLD + color.DARKCYAN + "Real Time:" + color.END + " %s secs"%"{:.2e}".format(t))
time_full = (time.time() - start_time)
print(color.BOLD + color.DARKCYAN + "Simulation Time:" + color.END + " %.7f secs"%(int(offset_full)+time_full))
print()
print(" ========= Number of reaction occurences ========= ")
for cat in range(len(Category)):
	print(Category[cat]+": %d"%counter_cat[cat])
