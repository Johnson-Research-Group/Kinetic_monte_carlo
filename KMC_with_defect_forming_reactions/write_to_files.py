import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

plt.figure(1)
plt.xlabel("Timesteps")

sizes = [8]
no_of_layers = 3
end_timestep=[16500]
no_of_C_atoms = [4608]


for row_ind in range(len(sizes)):
	print(" ")
	print("This is for size %d x %d"%(sizes[row_ind],sizes[row_ind]))
	print(" ")

	size=sizes[row_ind]
	
	f1=open("coordinates_%d_%d.cfg"%(size,no_of_layers),"r")
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
	f1=open("coordinates_O_%d_%d.cfg"%(size,no_of_layers),"r")
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
	coordinates_O_1 = np.stack( sorted(coordinates_O_1, key=lambda x: -x[2]), axis=0 )
	coordinates_O_3 = np.copy(coordinates)

	sites_O_1 = len(coordinates_O_1) #Epoxy
	sites_O_2 = len(coordinates_O_2) #Adatom
	sites_O_3 = len(coordinates_O_3) #In-plane

	Cat1=["Adsorption"]
	Cat2=["Thermal Desoprtion O Epoxy","O2 Recombination & Desorption", "Epoxy diffusion","Epoxy diffusion between layers", "1-Ep to 2-Ep","2-Ep to 3-Ep","1-Ep to 4-Ep", "1-Ep to Epoxy-Ether","Epoxy-Ether to 1-Ep","Epoxy-Ether to Epoxy-Ether-Epoxy (lactone-ether formation)"]
	Cat3=["Lactone-Ether CO formation","Ether-latone-ether CO formation","lactone-ether CO2 formation","ether-lactone-ether CO2 formation","lactone CO formation at VACANCY","lactone-epoxy CO formation at VACANCY","lactone CO2 formation at VACANCY","lactone-epoxy CO2 formation at VACANCY","Lactone-Ether to Ether-Lactone-Ether"]
	Cat4_=["O2 ZigZag Adsorption","O2 Armchair Adsorption","O ZigZag Adsorption","O ArmChair Adsorption","O Carbonyl Adsorption","O2 Split Diffusion ZZ1","O2 Split Diffusion AC1", "O2 Split Diffusion AC2","O2 Split Diffusion ZZ2", "O2 desorption ZZ", "O2 desorption AC", "No Epoxy ZZ CO formation", "1 Epoxy ZZ CO formation", "2 Epoxy ZZ CO formation", "1 carbonyl AC CO formation", "2 carbonyl AC CO formation", "1 carbonyl Epoxy AC CO formation", "2 carbonyl Epoxy AC CO formation", "1 carbonyl DZZ CO formation", "2 carbonyl DZZ CO formation", "3 carbonyl DZZ CO formation", "SBC CO formation", "SBC Epoxy CO formation"] 
	add_Cat4=[]

	for num in range(2,36):
		add_Cat4.append("%s Island Removal"%str(num)) 
	Cat4= Cat4_+add_Cat4
	Category1 = Cat1+Cat2
	Category2 = Cat3+Cat4
	Category = Category1+Category2

	time_=[]
	plot=[]
	range_ = range(0,end_timestep[row_ind],500)

	print("Timestep | Time Elapsed")
	for time_cat in range_:
		
		output = pd.read_pickle("Dataframe_%s_%d_%d.pkl"%(str(size), no_of_layers, time_cat)) #"All_flags_%d.pkl"
		flag, flag_O_1,flag_O_2, flag_O_3,offset_cat,offset_t, offset_step, offset_full, changing_av_ad, av_ad,changing_av_ep, av_ep, site_O_before, site_before, site_ad, site_ep, edge_C, cat = ((output.tail(1)).values.tolist())[0]

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
		color_O_3=[no_of_layers+5 if flag_O_3[i] == 'cyc-eth' else no_of_layers+7 for i in inds_to_keep_O_3]

		f = open("Sim_%s_%d_%d"%(str(size), no_of_layers, time_cat), "w")
		f.write("ITEM: TIMESTEP\n")
		f.write("%d\n"%time_cat)
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
		print("  %d      %.10f"%(time_cat,offset_t))
		time_.append(offset_t)
		plot.append(float(flag.count(['dg'])/no_of_C_atoms[row_ind]))

	plt.plot(time_, plot, label = "For size: %d"%size)

plt.xlabel("Time [seconds]")
plt.ylabel("Fraction of Carbon atoms removed")
plt.title("Rate of Carbon atom removal for %d layers"%no_of_layers)
plt.grid()
plt.legend()
plt.show()
