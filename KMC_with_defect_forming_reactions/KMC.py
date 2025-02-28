import numpy as np
import random
import math
import time
import itertools
from itertools import chain
import pandas as pd
import cantera as ct

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
	    return 1
	else:
		return 0

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

def armchair_pair(site):
	return next(p for p in near_neigh[site] if 'dg' not in flag[p] and [flag[k] for k in near_neigh[p]].count(['dg']) == 1)

def remove(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)
	for i in near_neigh[site]:
		Epoxy_Ad_C_delete(i)
		Lower_layer(i)

pos = [2.45953,2.12997,3.25363,3.68925]
def epoxy_curr_orient(site):
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

	val = [0]*5
	value=periodic(switch_site,site,coordinates_O_1,max_O_x+1.065,max_O_y+0.6149)
	if coordinates_O_1[switch_site][3] == coordinates_O_1[site][3]:
		val[1:4] = [isclose(value,pos[x], 0.1) for x in range(1,len(pos))]

		if isclose(value, 2.45953, 0.1) == 1 and len([x for x in epoxy_O_move[site] if x in epoxy_O_move[switch_site]])>0: val[0] = 1

	if val == [0]*5:
		val[4] = 1

	return val


def lactone_ether_group(site_C,flags):
	
	common_C = next(p for p in near_neigh[site_C]+[site_C] if 'ad' in flag[p] and flags in flag[p])

	stray_side_C = next(p for p in near_neigh[common_C] if 'stray' in flag[p])
	ether_C = next(p for p in near_neigh[common_C] if 'ether-side' in flag[p])
	lact_side_C = next(p for p in near_neigh[common_C] if 'lact-side' in flag[p])
	ether_C_side_C = next(p for p in near_neigh[ether_C] if p!=common_C and isclose(periodic(p,stray_side_C,coordinates,max_C_x+1.42,max_C_y+1.22975), 2.84, 0.1) == 1 and 'ether-side-side' in flag[p])
	cyc_eth_C = next(p for p in near_neigh[ether_C] if 'cyc-eth-C' in flag[p])
	return common_C, lact_side_C, ether_C, stray_side_C, ether_C_side_C, cyc_eth_C

def Epoxy_Ad_C_delete(C_atom):
	lay = int(coordinates[C_atom][3])
	flag_layers[lay]+=1

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
	cyc_eth_O = next(p for p in near_neigh[chosen] if flag_O_3[p] == 'cyc-eth')
	carb_cyc = [x for x in near_neigh[cyc_eth_O] if 'cyc-eth' in flag[x]]

	flag_O_3[cyc_eth_O] = 'na'
	cyc_ethers.remove(cyc_eth_O)
	for x in carb_cyc: flag[x] = ['carb_O']; flag_O_2[x] = 'carb_O'

	epoxy_remove_C(pivot_C)

def epoxy_diffusion_ep_2(site):
	common_C,lact_side_C,ether_C,stray_side_C,ether_C_side_C,cyc_eth_C  = lactone_ether_group(chosen, 'la-eth')
	list_ether_C = epoxy_list[ether_C]
	ether_O = [x[1] for x in list_ether_C if x[0]==common_C][0]
	ether_side_O=[x[1] for x in list_ether_C if x[0]==ether_C_side_C][0]

	flag[common_C],flag[lact_side_C],flag[ether_C],flag[stray_side_C],flag[ether_C_side_C], flag[cyc_eth_C],flag_O_2[common_C],flag_O_1[ether_O],flag_O_1[ether_side_O] = [['eth-la-eth','ad'], ['eth-la-eth','lact-side'], ['eth-la-eth','ether-side'], ['eth-la-eth','stray'], ['eth-la-eth','ether-side-side'],['eth-la-eth','cyc-eth-C'],'eth-la-eth','eth-la-eth','eth-la-eth']

	if flag[pivot_C]!= ['eth-la-eth','ether-side-side']:
			epoxy_remove_C(pivot_C)

def epoxy_diffusion_ep_3(site):
	flag_O_2[pivot_C]='carb_O'
	flag[pivot_C].remove('ep')
	flag[pivot_C].append('carb_O')

def epoxy_diffusion_ep_fast(site):
	switch_O_ep = [p for p in epoxy_O_move[site] if chosen in epoxy_O[p]][0]
	flag_O_1[site] = 'pr'
	flag_O_1[switch_O_ep] = 'ep'
	for i in epoxy_O[site]: flag[i] = ['pr']
	for i in epoxy_O[switch_O_ep]: flag[i] = ['ep']

def split_O2_adsorption(site):
	chosen_O = np.random.choice([x for x in epoxy_O_split[site] if (flag_O_1[x]=='pr' and [flag[epoxy_O[x][0]],flag[epoxy_O[x][1]]]==[['pr'],['pr']])])

	for i in epoxy_O[site]: epoxy_make_ep_C(i)
	flag_O_1[site] = 'ep'
	flag_O_1[chosen_O] = 'ep'

	for i in epoxy_O[chosen_O]: epoxy_make_ep_C(i)

def adsorption_ep(site):
	for i in epoxy_O[site]: epoxy_make_ep_C(i)
	flag_O_1[site] = 'ep'

def Recomb_O2_Desorp_ep(site):
	switch = int(np.random.choice((epoxy_curr_orient(site)[0])[1]))

	for j in [site, switch]:
		flag_O_1[j] = 'pr' 
		for i in epoxy_O[j]: epoxy_remove_C(i)

def adsorption_ad(site):
	flag_O_2[site] = 'ad'
	flag[site] = ['ad']

def epoxy_remove_C(site_C):
	if flag[site_C] == ['ep']:
			flag[site_C] = ['pr']
	else:
		flag[site_C].remove('ep')

def epoxy_make_ep_C(site_C):
	if flag[site_C] == ['pr']:
			flag[site_C] = ['ep']
	else:
		flag[site_C].append('ep')

def epoxy_move(site, switch_O):
	flag_O_1[site] = 'pr'
	for i in epoxy_O[site]:
		epoxy_remove_C(i)
	flag_O_1[switch_O]='ep'
	for i in epoxy_O[switch_O]:
		epoxy_make_ep_C(i)

def ep_diff_1_2(site):
	switch_O = np.random.choice(switch_O_1)
	epoxy_move(site, switch_O)

def ep_diff_2_3(site):
	switch_O = np.random.choice(switch_O_2)
	epoxy_move(site, switch_O)

def ep_diff_1_4(site):
	switch_O = np.random.choice(switch_O_3)
	epoxy_move(site, switch_O)

def pos1_to_Epoxy_Ether(site):
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
	ether_C = next(p for p in near_neigh[site] if 'eth' in flag[p] and 'la-vac' in flag[p])
	other_C = next(p for p in near_neigh[site] if p!=ether_C and 'dg' not in flag[p])
	ether_O = [x[1] for x in epoxy_list[site] if x[0]==ether_C][0]

	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

	flag[ether_C] =['pr']

def la_vac_epoxy_CO2(site): #lactone at vacancy to cycliic ether (CO2)
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
	flag[site] = ['carb_O2']
	flag_O_2[site] = 'carb_O2'

def O_edge_ad(site):
	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'

def O2_split_diff_ZZ1(site):
	chosen_C_site = next(i for i in ZigZag_neigh[site] if flag[i]==['pr'] and what_edge_carbon_is_it(i)[0]=="ZZ")
	flag[chosen_C_site] = ['carb_O']
	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[chosen_C_site] = 'carb_O'

def O2_split_diff_AC1(site):
	other_arm_chair = armchair_pair(site)
	flag[site] = ['carb_O']
	flag[other_arm_chair] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[other_arm_chair] = 'carb_O'

def O2_split_diff_AC2(site):
	other_arm_chair = next(i for i in Arm_chair_neigh[site] if flag[i]==['pr'] and what_edge_carbon_is_it(i)[0]=="AC")
	flag[site] = ['carb_O']
	flag[other_arm_chair] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	flag_O_2[other_arm_chair] = 'carb_O'

def O2_split_diff_ZZ2(site):
	chosen_C = next(i for i in near_neigh[site] if flag[i]==['pr'] and any(flag[k]==['pr'] and edge_C[k] !='edg' for k in near_neigh[i])==True)
	ep_C2 = next(i for i in near_neigh[chosen_C] if flag[i]==['pr'] and edge_C[i] !='edg')

	ep_O = [x[1] for x in epoxy_list[ep_C2] if x[0]==chosen_C][0]
	flag_O_1[ep_O] ='ep'
	flag[site] = ['carb_O']
	flag_O_2[site] = 'carb_O'
	for i in epoxy_O[ep_O]: epoxy_make_ep_C(i)

def O2_des(site):
	flag[site] =['pr']
	flag_O_2[site] = 'pr'

def no_ep_ZZ_CO(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def one_ep_ZZ_CO(site):
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
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def two_AC_CO(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def one_AC_ep_CO(site):
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
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def two_DZZ(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def three_DZZ(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def SBC_CO(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def SBC_ep_CO(site):
	Epoxy_Ad_C_delete(site)
	Lower_layer(site)

def add_upper_apron_x(Upper_apron_x, Upper_apron_y, Lower_apron_x, Lower_apron_y, coordinates_per):
	for l in layers:
		for i in Upper_apron_x[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0]-(max_C_x+0.71),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Upper_apron_y[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0],coordinates[i][1]-(max_C_y+1.22976),coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Lower_apron_x[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0]+(max_C_x+0.71),coordinates[i][1],coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
		for i in Lower_apron_y[l]:
			coordinates_per=coordinates_per+[[coordinates[i][0],coordinates[i][1]+(max_C_y+1.22976),coordinates[i][2],coordinates[i][3], coordinates[i][4]]]
	return coordinates_per

# parameters

size = 20
no_of_layers = 5

no_of_epoxies = 0

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
max_C_z = max(coordC_z)[0]

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
layernno_O = np.array([no_of_layers-1-int(math.ceil((x-1.29403)/3.4)) for x in coord_z]) 
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

#Constants
k_b_ev = 8.6173 * 10**(-5) #eV K^-1
k_b = 1.3806*(10**(-23)) #J/K
h = 6.62607*(10**(-34)) #Js
R_mol =8.314 #J/(Kmol) #m^3 Pa /(Kmol)
R = 259.84 #J/kgK
Av = 6.022*(10**23)

#Parameters
Surf_area = ((max_C_y*max_C_x)*(10**(-20))) #m^2
T_s = T_g = 1500
m_O = 2.657*(10**(-26)) #kg
m_C = 1.99 * 10**(-26) #kg
M_O = 16*(10**-3) #kg/mol
rho_o = 1.293 #kg/m^3
rho = np.exp(0)*rho_o

Altitude = 30.8
Mach = 6.31
vel = np.sqrt((3*k_b*T_g)/m_O) #m/s
delta_t = 10**(-10) #add this time for adsorption rate

stick_coeff = 0.85 #for 5eV, Increases with decreasing beam energy

# J --> kg⋅m2⋅s−2
# Pa --> kg⋅m−1⋅s−2
P_stag = 10000
gas1 = ct.Solution('gri30.yaml')
gas1.TPX = T_g, P_stag, 'O2:0.21, N2:0.78, Ar:0.09'
gas1.equilibrate('TP')
X = gas1['O','O2'].X
P_O = X[0]*P_stag
P_O2 = X[1]*P_stag

print("For size %d"%size)
print("For number of layers %d"%no_of_layers)
print("For temperature: %d"%T_s)
print("For Pressure: %d"%P_stag)
print("Partial Pressure O: %s"%str(P_O))
print("Partial Pressure O2: %s"%str(P_O2))

number_density = P_O/(k_b*T_g)

Cat1=["Adsorption","Split O2 Adsorpiton"]
Cat2=["O2 Recombination & Desorption", "Epoxy diffusion Fast","Epoxy diffusion between layers Fast","Epoxy diffusion","Epoxy diffusion between layers", "1-Ep to 2-Ep","2-Ep to 3-Ep","1-Ep to 4-Ep", "1-Ep to Epoxy-Ether","Epoxy-Ether to 1-Ep","Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation)"]
Cat3=["Lactone-Ether CO formation","Ether-latone-ether CO formation","lactone-ether CO2 formation","ether-lactone-ether CO2 formation","lactone CO formation at VACANCY","lactone-epoxy CO formation at VACANCY","lactone CO2 formation at VACANCY","lactone-epoxy CO2 formation at VACANCY","Lactone-Ether to Ether-Lactone-Ether"]
Cat4_=["O2 ZigZag Adsorption","O2 Armchair Adsorption","O ZigZag Adsorption","O ArmChair Adsorption","O Carbonyl Adsorption","O2 Split Diffusion ZZ1","O2 Split Diffusion AC1", "O2 Split Diffusion AC2","O2 Split Diffusion ZZ2", "O2 desorption ZZ", "O2 desorption AC", "No Epoxy ZZ CO formation", "1 Epoxy ZZ CO formation", "2 Epoxy ZZ CO formation", "1 carbonyl AC CO formation", "2 carbonyl AC CO formation", "1 carbonyl Epoxy AC CO formation", "2 carbonyl Epoxy AC CO formation", "1 carbonyl DZZ CO formation", "2 carbonyl DZZ CO formation", "3 carbonyl DZZ CO formation", "SBC CO formation", "SBC Epoxy CO formation"] 
add_Cat4=[]
for num in range(2,36):
	add_Cat4.append("%s Island Removal"%str(num)) 
Cat4= Cat4_+add_Cat4
Category1 = Cat1+Cat2
Category2 = Cat3+Cat4
Category = Category1+Category2

N1 = len(Category1)
N2 = len(Category2)

k1 = np.empty(N1, dtype=object) 
k2 = np.empty(N2, dtype=object) 

#========== Surface ===========

#GAS REACTIONS

A_eff = 2.6199 * (10**-20)

#----- Adsorption
k1[0] = ((P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)) #*k_inc #Adsorption #1/s #-------> Microscopic modeling and optimal operation of plasma enhanced atomic layer deposition (Eq 3)
k1[1] = ((P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g))*np.exp(-1.7/(k_b_ev*T_s))  #O2 split adsorption # 1/s # change to O2 values

#----- Desorption
#k1[2]=k1[3] = 0 #(k_b*T_s/h)*np.exp(-1.91/(k_b_ev*T_s)) #Thermal Desorption Epoxy #1/s 
k1[2] = 2.029*(10**(-17))*number_density*np.exp(-1.130017/(k_b_ev*T_s)) #2-Ep to O2 Recombination & Desorption #1/s

#DIFFUSION

#------ Epoxy diffusion
k1[3]=k1[5] = (10**13)*np.exp(-0.737/(k_b_ev*T_s)) #Epoxy Diffusion #1/s #-------> Dual Path Mechanism in the Thermal Reduction of Graphene Oxide, Stefano Fabris, Tao Sun
k1[4] =k1[6]= (10**13)*np.exp(-0.737/(k_b_ev*T_s)) #Epoxy Diffusion between layers #0.8, 0.9
k1[7] = (10**13)*np.exp(-0.66/(k_b_ev*T_s)) # 1-Ep to 2-Ep #1/s
k1[8] = (10**13)*np.exp(-1.13/(k_b_ev*T_s)) # 2-Ep to 3-Ep #1/s
k1[9] = (10**13)*np.exp(-1.34/(k_b_ev*T_s)) # 1-Ep to 4-Ep #1/s
k1[10] = (10**13)*np.exp(-1.34/(k_b_ev*T_s)) # 1-Ep to Epoxy-Ether #1/s #-1.34
k1[11] = (10**13)*np.exp(-0.95/(k_b_ev*T_s)) # Epoxy-Ether to 1-Ep #1/s #-0.95
k1[12] = (10**13)*np.exp(-1.1/(k_b_ev*T_s)) # Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation) #1/s

#SURFACE REACTIONS
#------ CO
#Below 3ev
k2[0] = (k_b*T_s/h)*np.exp(-0.97/(k_b_ev*T_s)) #Lactone-Ether CO formation
k2[1] = (k_b*T_s/h)*np.exp(-0.5/(k_b_ev*T_s)) #Ether-latone-ether CO formation
#k2[2] = (k_b*T_s/h)*np.exp(-2.6/(k_b_ev*T_s)) #lactone CO formation
k2[4] = (k_b*T_s/h)*np.exp(-2.98/(k_b_ev*T_s)) #lactone CO formation at VACANCY #2.98
k2[5] = (k_b*T_s/h)*np.exp(-1.93/(k_b_ev*T_s)) #lactone-epoxy CO formation at VACANCY #1.93

#------ CO2
#Below 3eV
#k3[5] = (k_b*T_s/h)*np.exp(-1.84/(k_b_ev*T_s)) #lactone CO2 formation
k2[2] = (k_b*T_s/h)*np.exp(-0.61/(k_b_ev*T_s)) #lactone-ether CO2 formation
k2[3] = (k_b*T_s/h)*np.exp(-0.6/(k_b_ev*T_s)) #ether-lactone-ether CO2 formation 
k2[6] = (k_b*T_s/h)*np.exp(-2.46/(k_b_ev*T_s)) #lactone CO2 formation at VACANCY #2.46
k2[7] = (k_b*T_s/h)*np.exp(-0.58/(k_b_ev*T_s)) #lactone-epoxy CO2 at VACANCY #0.58
k2[8] = (k_b*T_s/h) #Lactone-Ether + epoxy to Ether-Lactone-Ether
#============= Edge ===============

#GAS REACTIONS
k2[9] = (P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g)  # O2 Adsorption ZigZag #1/s 
k2[10] = (P_O2*A_eff)/np.sqrt(2*np.pi*2*m_O*k_b*T_g)*np.exp(-0.311/(k_b_ev*T_s)) # O2 Adsorption ArmChair #1/s
k2[11] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)  # O ZigZag Adsorption
k2[12] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g) # O ArmChair Adsorption
k2[13] = (P_O*A_eff)/np.sqrt(2*np.pi*m_O*k_b*T_g)  # O Carbonyl Adsorption
k2[14] = 3.6*(10**12)*np.exp(-0.041/(k_b_ev*T_s)) #O2 Split Diffusion ZZ1
k2[15] = 1.5*(10**12)*np.exp(-0.694/(k_b_ev*T_s)) #O2 Split Diffusion AC1
k2[16] = 3.5*(10**12)*np.exp(-0.238/(k_b_ev*T_s)) #O2 Split Diffusion AC2
k2[17] = 1.2*(10**12)*np.exp(-1.513/(k_b_ev*T_s)) #O2 Split Diffusion ZZ2
k2[18] = 1.3*(10**14)*np.exp(-1.42/(k_b_ev*T_s)) #O2 desorption ZZ
k2[19] = 6.4*(10**14)*np.exp(-1.104/(k_b_ev*T_s)) #O2 desorption AC
k2[20] = 1.2*(10**16)*np.exp(-3.622/(k_b_ev*T_s)) #No Epoxy ZZ CO formation
k2[21] = 1.2*(10**16)*np.exp(-2.3/(k_b_ev*T_s)) #1 Epoxy ZZ CO formation
k2[22] = 1.2*(10**16)*np.exp(-2.3/(k_b_ev*T_s)) #2 Epoxy ZZ CO formation
k2[23] = 1.0*(10**13)*np.exp(-2.743/(k_b_ev*T_s)) #1 carbonyl AC CO formation
k2[24] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) #2 carbonyl AC CO formation
k2[25] = 1.0*(10**13)*np.exp(-1.829/(k_b_ev*T_s)) #1 carbonyl Epoxy AC CO formation
k2[26] = 1.0*(10**13)*np.exp(-1.136/(k_b_ev*T_s)) #2 carbonyl Epoxy AC CO formation
k2[27] = 1.0*(10**13)*np.exp(-2.743/(k_b_ev*T_s)) #1 carbonyl DZZ CO formation
k2[28] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) #2 carbonyl DZZ CO formation
k2[29] = 1.0*(10**13)*np.exp(-1.704/(k_b_ev*T_s)) #3 carbonyl DZZ CO formation
k2[30] = 1.0*(10**13)*np.exp(-2.17/(k_b_ev*T_s)) #SBC CO formation
k2[31] = 1.0*(10**13)*np.exp(-0.447/(k_b_ev*T_s)) #SBC Epoxy CO formation

#Island Removal
for num in range(2,36):
	k2[num+31-1] = 1.0*(10**13)*np.exp(-(0.0751*num)/(k_b_ev*T_s)) # num carbon island

#==========================================================

k1 = k1.astype('float')
k2 = k2.astype('float')
k_const = np.concatenate(((k1,k2)), axis=0)

init_time1 = time.time()

cov_C = ['uncov' if coordinates[i][3]==0 else 'cov' for i in range(sites)]
cov_ep = ['uncov' if coordinates_O_1[i][3]==0 else 'cov' for i in range(sites_O_1)]
edge_C = ['pr' for i in range(sites)]
offset_cat=0
offset_t=0
offset_step=0
offset_full=0
site_O_before='na'
site_before='na'
switch_ad_before = [0]*no_of_layers
checker=0

changing_av_ad = np.zeros((sites_O_2,N2))
changing_av_ep = np.zeros((sites_O_1,N1))

sites_av_ad = [[] for i in range(N2)]
sites_av_ep = [[] for i in range(N1)]

flag_layers = [0 for i in range(no_of_layers)]
flag_check = [0 for i in range(no_of_layers)]

top_layer_O_1 = [i for i in range(sites) if coordinates_O_1[i][3]==0]

for i in top_layer_O_1:
	changing_av_ep[i][0]=1
	changing_av_ep[i][1]=1

flag = [['pr'] for i in range(sites)]
flag_O_1 = ['pr' for i in range(sites_O_1)]
flag_O_2 = ['pr' for i in range(sites_O_2)]
flag_O_3 = ['na' for i in range(sites_O_3)]

sur_av_ep = np.array([len((np.where(changing_av_ep[:,n]==1)[0]).tolist()) for n in range(N1)])
sur_av_ad = np.array([len((np.where(changing_av_ad[:,n]==1)[0]).tolist()) for n in range(N2)])
sur_av = sur_av_ep.tolist()+sur_av_ad.tolist()

counter_cat = [float(0)]*len(sur_av)

Elapsed_time =[]

init_time2 = time.time()

C_remove_end1 = 0
C_remove_end2 = 0
cat=0

initialize_time = init_time2-init_time1

Epoxy_time_list = []
Adatom_time_list = []
compile_av_event_list = []
choosing_ev_list = []
choosing_site_list = []
perform_event_list = []
saving_file_list = []

def epoxy_rit(site_O):
	before_ep_arr = (np.copy(changing_av_ep[site_O]))
	changing_av_ep[site_O]=np.array(n1)
	diff = (changing_av_ep[site_O]-before_ep_arr)
	sur_av_ep[:] = sur_av_ep[:] + diff

def adatom_rit(site):
	before_ad_arr = np.copy(changing_av_ad[site])
	changing_av_ad[site]=np.array(n2)
	diff = changing_av_ad[site]-before_ad_arr
	sur_av_ad[:] = sur_av_ad[:] + diff

cyc_ethers = []

#For re-running code from previous run

# output3 = pd.read_pickle("Df_temp%d_19_7_22_%dPa_%dK_%s_%d_%d.pkl"%(set_,P_stag,T_g,str(size), no_of_layers, restart_timestep))
# flag, flag_O_1, flag_O_2, flag_O_3, offset_cat,offset_t, offset_step, offset_full, changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat,flag_layers = ((output3.tail(1)).values.tolist())[0]

# flag_check = [1 if flag_layers[lay]>int(no_C_each_layer/2) else 0 for lay in range(no_of_layers)]

# time
t = offset_t

#KMC
time_step = int(offset_step)
time_cat = int(offset_cat)

while cat!=-1: #flag.count('dg') != sites 

	time1 = time.time()

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

	[ind_epoxy_pr.append(x) if (flag_O_1[x] == 'pr' and [flag[epoxy_O[x][0]],flag[epoxy_O[x][1]]] ==[['pr'],['pr']]) and cov_ep[x]=='uncov' else ind_epoxy_ep.append(x) if flag_O_1[x] == 'ep' else ind_epoxy_eth_ep.append(x) if (flag_O_1[x] == 'ep-eth' and any('epoxy' in flag[i] for i in epoxy_O[x])==True) else ind_epoxy_set.append(x) for x in ind_epoxy]

	n1 = [0]*N1
	[epoxy_rit(site_O) for site_O in ind_epoxy_set]

	#time0_2 = time.time()

	ind_epoxy_pr_nor = []
	ind_epoxy_pr_split = []
	[ind_epoxy_pr_split.append(x) if any((flag_O_1[p] == 'pr' and [flag[epoxy_O[p][0]],flag[epoxy_O[p][1]]] ==[['pr'],['pr']]) and cov_ep[p]=='uncov' for p in epoxy_O_split[x])==True else ind_epoxy_pr_nor.append(x) for x in ind_epoxy_pr]

	#time0_3 = time.time()

	n1[0]=1
	[epoxy_rit(site_O) for site_O in ind_epoxy_pr_nor]
	n1[1]=1
	[epoxy_rit(site_O) for site_O in ind_epoxy_pr_split]
	
	for site_O in ind_epoxy_ep:
		n1 = [0]*N1
		curr_list, curr_orient = epoxy_curr_orient(site_O)
		covered_val = cov_ep[site_O]

		n1[2]= 1 if curr_orient[1]==1 else 0 #2-Ep to O2 Recombination & Desorption

		C_move=[-1]*6
		if all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and all(flag[x]==['pr'] for x in near_neigh_8[site_O])==True:
			n1[3] = 1
			#n1[4] = 1 if any(cov_C[x]=='cov' for x in epoxy_C_move[site_O])==True else 0

		elif all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and all('ep' not in flag[x] for x in near_neigh_8[site_O])==True and all(flag[x]==['pr'] for x in near_neigh_8[site_O])!=True:
			occupied = [x for x in epoxy_C_move[site_O] if any(i in ['eth-la-eth','carb_O','carb_O2','la-vac','ep-eth','la-eth'] for i in flag[x])==True]
			n1[5] = 1 if len(occupied)<4 else 0
			#n1[6] = 1 if blocked<4 and any(cov_C[x]=='cov' for x in possible)==True else 0

		elif all(flag[p]==['ep'] for p in epoxy_O[site_O])==True and any('ep' in flag[x] for x in near_neigh_8[site_O])==True:
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
			n1[5] = 1 if (len(C_move)+len(occupied))<4 else 0
			#n1[6] = 1 if (len(C_move)+len(occupied))<4 and any(cov_C[x]=='cov' for x in possible)==True else 0

			n1[7] = 1 if a==1 else 0 # 1-Ep to 2-Ep 
			n1[8] = 1 if b==1 else 0 # 2-Ep to 3-Ep 
			n1[9] = 1 if c==1 else 0 # 1-Ep to 4-Ep 
			n1[10] = 1 if d==1 else 0 # 1-Ep to Epoxy-Ether

		else:
			pass

		epoxy_rit(site_O)

		if all(flag[p] == ['dg'] for p in epoxy_C_move[site_O])==True:
			for i in epoxy_O[site_O]:
				Epoxy_Ad_C_delete(i)
				Lower_layer(i)

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

	#SURFACE REACTIONS
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

	for site in ind_ad_edg:
		n2 = [0]*N2
		check_edge_C = 0
		#==============================================================================#
		edge_C_site_type,edge_C_list = what_edge_carbon_is_it(site)

		#Islands
		#print("================ Islands ==================")
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

		#Rate constants
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

			# Diffusion
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
					#n2[14] = 1 if any(flag[i]==['pr'] and any(p in ["ZZ", "AC"] for p in what_edge_carbon_is_it(i)[0])==True for i in ZigZag_neigh[site])==True else 0
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

	sur_av = sur_av_ep.tolist()+sur_av_ad.tolist()

	rates = np.absolute(np.einsum('i,i->i',k_const,sur_av))
	R = np.einsum('i,i->',k_const,sur_av)

	if R != 0:
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


	if cat<N1 and cat>1: # (Greater than 0 and not -1 because no prep needed for adsorption)
		site_O = np.random.choice(np.where(changing_av_ep[:,cat]==1)[0])
		site_O_before=site_O

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

	elif cat>(N1-1) and cat< len(sur_av)-1:
		site = np.random.choice(np.where(changing_av_ad[:,cat-N1]==1)[0])
		site_before=site
		#print("This is the site: %d"%site)

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

			# Diffusion
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

	elif cat == 7: #1-Ep to 2-Ep
		ep_diff_1_2(site_O)

	elif cat == 8:#2-Ep to 3-Ep
		ep_diff_2_3(site_O)

	elif cat == 9: #1-Ep to 4-Ep
		ep_diff_1_4(site_O)

	elif cat == 10: #1-Ep to Epoxy-Ether
		pos1_to_Epoxy_Ether(site_O)

	elif cat == 11: #Epoxy-Ether to 1-Ep
		Epoxy_Ether_to_pos1(site_O)

	elif cat == 12: #Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation)
		Ep_Eth__Ep_Eth_Ep(site_O)

	elif cat == 13: #Lactone-Ether CO formation
		la_eth_CO(site)

	elif cat == 14: #Ether-latone-ether CO formation
		eth_la_eth_CO(site)

	elif cat ==15: #lactone-ether CO2 formation
		la_eth_CO2(site)

	elif cat == 16:#ether-lactone-ether CO2 formation
		eth_la_eth_CO2(site)

	elif cat == 17: #lactone CO formation at VACANCY
		la_vac_CO(site)

	elif cat == 18: #lactone-epoxy CO formation at VACANCY
		la_vac_epoxy_CO(site)

	elif cat == 19: #lactone CO2 formation at VACANCY
		la_vac_CO2(site)

	elif cat == 20: #lactone-epoxy CO2 formation at VACANCY
		la_vac_epoxy_CO2(site)

	elif cat == 21:
		la_eth_to_eth_la_eth(site)

	if cat == 22 or cat == 23: #O2 ZigZag & Armchair Adsorption 
		O2_edge_ad(site)

	elif cat ==24 or cat == 25 or cat == 26: #O ZigZag Adsorption & O ArmChair Adsorption & O Carbonyl Adsorption
		O_edge_ad(site)

	elif cat == 27: #O2 Split Diffusion ZZ1
		O2_split_diff_ZZ1(site)

	elif cat == 28: #O2 Split Diffusion AC1
		O2_split_diff_AC1(site)

	elif cat == 29: #O2 Split Diffusion AC2
		O2_split_diff_AC2(site)

	elif cat == 30: #O2 Split Diffusion ZZ
		O2_split_diff_ZZ2(site)
		
	elif cat == 31: #O2 Desorption ZZ
		O2_des(site)

	elif cat==32: #O2 Desorption AC
		O2_des(site)

	elif cat==33: # No Epoxy ZZ CO formation
		no_ep_ZZ_CO(site)

	elif cat==34: # 1 Epoxy ZZ CO formation
		one_ep_ZZ_CO(site)

	elif cat==35: # 2 Epoxy ZZ CO formation
		two_ep_ZZ_CO(site)

	elif cat==36: # 1 carbonyl AC CO formation
		one_AC_CO(site)

	elif cat==37: # 2 carbonyl AC CO formation
		two_AC_CO(site)

	elif cat==38: # 1 carbonyl Epoxy AC CO formation
		one_AC_ep_CO(site)

	elif cat==39: # 2 carbonyl Epoxy AC CO formation
		two_AC_ep_CO(site)

	elif cat==40: # 1 carbonyl DZZ CO formation
		one_DZZ(site)

	elif cat==41: # 2 carbonyl DZZ CO formation
		two_DZZ(site)

	elif cat==42: # 3 carbonyl DZZ CO formation
		three_DZZ(site)

	elif cat==43: # SBC CO formation
		SBC_CO(site)

	elif cat==44: # SBC Epoxy CO formation
		SBC_ep_CO(site)

	elif cat>=45 and cat<(36+45):
		for elem in chain:
			Epoxy_Ad_C_delete(elem)
			Lower_layer(elem)

	else: #ZigZag Scatter & Armchair Scatter & Carbonyl O Scatter
		pass

	#Cyclic Ethers

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

	if time_step%1000000==0:
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3, time_cat,t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])
		df1_4.to_pickle("./Df_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))
		print(" ========= KMC Iteration = %d ==========" %time_step)
	
	if cat==-1:
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3, time_cat,t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])#, columns =['Carbon Flags','Epoxy Flags''Carbonyl 
		df1_4.to_pickle("./Df_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))

	if (time.time() - start_time)/3600 >69 and checker==0:
		checker=1
		df1_4 = pd.DataFrame([[flag, flag_O_1, flag_O_2, flag_O_3, time_cat,t,time_step, int(offset_full)+(time.time() - start_time), changing_av_ad, changing_av_ep, site_O_before, site_before, sur_av_ad, sur_av_ep, edge_C, cat,prob,cov_C,cov_ep,counter_cat]])
		df1_4.to_pickle("./Df_restart_%dPa_%dK_%s_%d_%d.pkl"%(P_stag,T_g,str(size), no_of_layers, time_step))

	if cat not in [3,4,5,6,len(sur_av)-1]:
		time_cat +=1

	if R!=0:
		r1 = random.random()
		t += -(np.log(r1))/(R*1.0)
	else:
		pass
	
	#timestep_end = time.time()
	#print("Time for timestep: %s seconds"%str(timestep_end -time1))
	#print("Time Elapsed: %s mins"%str((time.time() - start_time)/60))
	#print(" ========= t_step = %d ==========" %time_step)

	time_step +=1


print("Number of Epoxies: %d" %flag_O_1.count('ep'))
print(color.BOLD + color.DARKCYAN + "Real Time:" + color.END + " %s secs"%"{:.2e}".format(t))
time_full = (time.time() - start_time)
print(color.BOLD + color.DARKCYAN + "Simulation Time:" + color.END + " %.7f secs"%(int(offset_full)+time_full))
