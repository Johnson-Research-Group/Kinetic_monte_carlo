from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import time
import csv
import random
import math

from ctypes import *

f=open("FloatingCO_15A_145eV_defectl4.txt","a+")
for s in range(145,146):
    for itr in range (1,21):
        for itr1 in range(1,2):
            from lammps import lammps, PyLammps
            lmp = lammps()
            L = PyLammps(ptr=lmp)
            L.atom_modify("map yes")
            L.read_restart("heated_300K_defectlayer4.txt")
            L.variable("length equal lx")
            L.variable("l equal ${length}/30")
            L.variable("zlength equal lz")
            L.variable("layers equal 5")
            L.variable("zlat equal (${zlength}-135)/(${layers}-1)+3.35")
            L.variable("s equal ${zlat}")
            L.variable("a1 equal 3*$l")
            L.variable("a2 equal (sqrt(3))*$l")
            L.variable("a3 equal $s*2")

        #Change height of oxygen for randomness
            
            #randOposition=round(random.uniform(1,11),3)
            #L.variable("ho equal 2*${a3}+7+%.3f"%(randOposition))
            #heightofO = float(float(7)+randOposition)
            dheightO=(float(itr)-10)/1000
            L.variable("ho equal 2*${a3}+15+(%.3f)"%(dheightO))
            heightofO = float(float(15)+dheightO)
            print(heightofO)

        #Change position of oxygen in square for randomness

            #randOposX=round(random.uniform(0,3.16),3)
            #randOposY=round(random.uniform(0,4.26),3)
            #L.variable("o1 equal (5*${a1})+2.5*$l+%.3f"%randOposX)
            #L.variable("o2 equal (5*${a2})+((sqrt(3))/4)*$l+%.3f"%randOposY)
            
            #a2=math.sqrt(3)*1.42
            #Xposo = float((5*4.26)+(2.5*1.42)+float(randOposX))
            #Yposo = float((5*a2)+((math.sqrt(3)*1.42)/4)+float(randOposY))
            
            if s==0:
                eV = 5
            else:
                eV = s*1

            print(eV)
            L.variable("T equal %d" % eV)
            L.file("in.graphene11")

        # Finding the maximum and minimum zcoord of atoms in each carbon layer
            L.run(0)
            l1id=[]
            l2id=[]
            l3id=[]
            l4id=[]
            l5id=[]
            zall=[]

            for i in range(0,3000):
                pz=L.atoms[i].position
                zall.append(pz[2])

            for i in range(0,3000):
                if zall[i]>11.725:
                    l1id.append(zall[i])
                elif 8.375<zall[i]<11.725:
                    l2id.append(zall[i])
                elif 5.025<zall[i]<8.375:
                    l3id.append(zall[i])
                elif 1.675<zall[i]<5.025:
                    l4id.append(zall[i])
                else:
                    l5id.append(zall[i])

            for r in range(1,6):
                exec "l%dmaxz=max(l%did)+0.7"%(r,r)
                exec "l%dminz=min(l%did)-0.7"%(r,r)

        #RUN CODE 

            L.run(2000)

            # NUMBER OF BONDS OF OXYGEN ATOM
            f1=open("bonds.reaxc", "r")    
            file1=f1.readlines()
            count=0
            for line in file1:
                count += 1
            linesno = count

            x1=file1[(count-3002):(count-1)]
            idnumber=[]
            numberofbonds=[]
            bid1=[]
            bid2=[]
            bid3=[]

            for lines in x1:
                idnumber.append(lines.split()[0])
                numberofbonds.append(lines.split()[2])
                bid1.append(lines.split()[3])
                bid2.append(lines.split()[4])
                bid3.append(lines.split()[5])
                

            f1.close()
            
            #a=idnumber.index("3001")
            a=idnumber.index("917")

            def maybe_float(s):
                try:
                    return float(s)
                except (ValueError, TypeError):
                    return s
            pno = [maybe_float(i) for i in idnumber]
            nob = [maybe_float(i) for i in numberofbonds]
            bid_1 = [maybe_float(i) for i in bid1]
            bid_2 = [maybe_float(i) for i in bid2]
            bid_3 = [maybe_float(i) for i in bid3]
            b=nob[a]
            print(b)

        # Z LOCATION OF OXYGEN ATOM
            #p=L.atoms[3000].position
            p=L.atoms[916].position
            zo=p[2]

        # ID of atom(s) attached to oxygen
            #print(bid_1)
            #c1id_to_oxygen = bid_1[a]
            #if c1id_to_oxygen==0:
            #    pass
            #else:
            #    k = pno.index(c1id_to_oxygen)
            #    bc1id_to_oxygen=nob[k]

            #if bid2!="0":
            #   c2id_to_oxygen = bid_2[a]
            #   if c2id_to_oxygen==0:
            #       pass
            #   else:
            #      k1 = pno.index(c2id_to_oxygen)
            #       bc2id_to_oxygen=nob[k1] 

            if b==0:
                pass
            elif b==1:
                c1id_to_oxygen = bid_1[a]
                k = pno.index(c1id_to_oxygen)
                bc1id_to_oxygen=nob[k]
           
            elif b==2:
                c1id_to_oxygen = bid_1[a]
                k = pno.index(c1id_to_oxygen)
                bc1id_to_oxygen=nob[k]
                c2id_to_oxygen = bid_2[a]
                k1 = pno.index(c2id_to_oxygen)
                bc2id_to_oxygen=nob[k1] 
            elif b==3:
                c1id_to_oxygen = bid_1[a]
                k = pno.index(c1id_to_oxygen)
                bc1id_to_oxygen=nob[k]
                c2id_to_oxygen = bid_2[a]
                k1 = pno.index(c2id_to_oxygen)
                bc2id_to_oxygen=nob[k1] 
                c3id_to_oxygen = bid_3[a]
                k2 = pno.index(c3id_to_oxygen)
                bc3id_to_oxygen=nob[k2] 
            else:
                pass

        # Layers ID AND Z LOCATION OF LOWEST CARBON ATOM IN ALL LAYER
            d={}
            for j in range(1,6):
                d["f"+str(j+1)]=open("Layer%d.cfg"%j, "r") 
                d["file"+str(j+1)]=d["f"+str(j+1)].readlines()
                d["count"+str(j+1)]=0
                for line in d["file"+str(j+1)]:
                    d["count"+str(j+1)] += 1
                d["linesno"+str(j+1)]=d["count"+str(j+1)]
                d["x"+str(j+1)]=d["file"+str(j+1)][(d["count"+str(j+1)]-600):(d["count"+str(j+1)])]
                
                d["layer"+str(j)+"id"] = []
                d["zcoordinate"+str(j)]=[]

                for lines in d["x"+str(j+1)]:
                    d["zcoordinate"+str(j)].append(lines.split()[2])
                    d["layer"+str(j)+"id"].append(lines.split()[0])

                d["zlist"+str(j)] = [maybe_float(i) for i in d["zcoordinate"+str(j)]]
                d["layer"+str(j)+"_id"] = [maybe_float(i) for i in d["layer"+str(j)+"id"]]
                d["zcmin"+str(j)] = float(min(d["zlist"+str(j)]))

        # CARBON ATOMS IN LAYER1 WITH 4 BONDS
          
            res = [] 
            i = 0
            while (i < len(d["layer1_id"])): 
                if (pno.count(d["layer1_id"][i]) > 0): 
                    res.append(i) 
                i += 1

            #print(res)
            
            cbonds = [nob[i] for i in res]
            cb1_4=cbonds.count(4)
            cb1_2=cbonds.count(2)
            cb1_1=cbonds.count(1)
                
                #print(len(res))
                #print(cb1_4)

            tag=[" "," "," "," "]
            print(zo)
            print(b)

        # CATEGORIES
            if zo>15:
                tag[0]="Oxygen scattered back" #Oxygen is scattered back
        #Oxygen just above layer 1
            elif l1maxz<zo<15:
                tag[0]= "Oxygen attached to layer 1"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                    if cb1_4==1 and cb1_2==0 and zbetl12=="no": 
                        tag[2] = "No broken layers" #Oxygen is attached to layer1 by one bond with no bonds broken in Layer 1
                    elif cb1_4==0 and l1minz<d["zcmin1"]<l2maxz and zbetl12=="yes":
                        tag[2] = "Break in layer 1" #Break in Layer 1
                    else:
                        pass

                elif b==2:
                    tag[1]="Epoxy group"
                    if cb1_4==2 and cb1_2==0:
                        tag[2]="No broken C-C bond in Epoxy" # Epoxy group 
                    elif cb1_4==0 and cb1_2==0:
                        tag[2]="broken C-C bond in Epoxy" # Epoxy group (broken C-C bond)
                    else:
                        pass
                else:
                    pass
            
            elif l1minz<zo<l1maxz:
                tag[0]="Oxygen attached to layer 1"
                if (b==3 or b==2):
                    tag[1]="Oxygen part of graphene sheet"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                else:
                    pass

            elif l2minz<zo<l2maxz:
                tag[0]="Oxygen attached to layer 2"
                if (b==3 or b==2):
                    tag[1]="Oxygen part of graphene sheet"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                else:
                    pass
         

            elif l3minz<zo<l3maxz:
                tag[0]="Oxygen attached to layer 3"
                if (b==3 or b==2):
                    tag[1]="Oxygen part of graphene sheet"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                else:
                    pass

            elif l4minz<zo<l4maxz:
                tag[0]="Oxygen attached to layer 4"
                if (b==3 or b==2):
                    tag[1]="Oxygen part of graphene sheet"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                else:
                    pass

            elif l5minz<zo<l5maxz:
                tag[0]="Oxygen attached to layer 5"
                if (b==3 or b==2):
                    tag[1]="Oxygen part of graphene sheet"
                if b==1:
                    tag[1]="Oxygen attached with single bond to layer"
                else:
                    pass

            elif l2maxz<zo<l1minz:
                if b==0:
                    tag[0]="Between layers 1&2"
                    tag[1]="Floating Oxygen atom"
                elif b==1 and bc1id_to_oxygen==1:
                    tag[0]="Between layers 1&2"
                    tag[1]="Floating C-O bond"
                elif b==1 and (c1id_to_oxygen in d["layer1_id"]):
                    tag[0]="Oxygen attached to layer 1"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==1 and (c1id_to_oxygen in d["layer2_id"]):
                    tag[0]="Oxygen attached to layer 2"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==2 and (c1id_to_oxygen in d["layer1_id"]) and (c2id_to_oxygen in d["layer1_id"]):
                    tag[0]="Oxygen attached to layer 1"
                    tag[1]="Epoxy group"
                elif b==2 and (c1id_to_oxygen in d["layer2_id"]) and (c2id_to_oxygen in d["layer2_id"]):
                    tag[0]="Oxygen attached to layer 2"
                    tag[1]="Epoxy group"
                elif b==2 and (((c1id_to_oxygen in d["layer1_id"]) and (c2id_to_oxygen in d["layer2_id"])) or ((c1id_to_oxygen in d["layer2_id"]) and (c2id_to_oxygen in d["layer1_id"]))):
                    tag[0]="Between layers 1&2"
                    tag[1]="Bridge Oxygen"
                elif b==3:
                    tag[0]="Oxygen attached to layer 2"
                    tag[1]="Epoxy group shifting betwen 2/3 bonds"
                else:
                    tag[0]="Between layers 1&2"

            elif l3maxz<zo<l2minz:
                if b==0:
                    tag[0]="Between layers 2&3"
                    tag[1]="Floating Oxygen atom"
                elif b==1 and bc1id_to_oxygen==1:
                    tag[0]="Between layers 2&3"
                    tag[1]="Floating C-O bond"
                elif b==1 and (c1id_to_oxygen in d["layer2_id"]):
                    tag[0]="Oxygen attached to layer 2"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==1 and (c1id_to_oxygen in d["layer3_id"]):
                    tag[0]="Oxygen attached to layer 3"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==2 and (c1id_to_oxygen in d["layer2_id"]) and (c2id_to_oxygen in d["layer2_id"]):
                    tag[0]="Oxygen attached to layer 2"
                    tag[1]="Epoxy group"
                elif b==2 and (c1id_to_oxygen in d["layer3_id"]) and (c2id_to_oxygen in d["layer3_id"]):
                    tag[0]="Oxygen attached to layer 3"
                    tag[1]="Epoxy group"
                elif b==2 and (((c1id_to_oxygen in d["layer2_id"]) and (c2id_to_oxygen in d["layer3_id"])) or ((c1id_to_oxygen in d["layer3_id"]) and (c2id_to_oxygen in d["layer2_id"]))):
                    tag[0]="Between layers 2&3"
                    tag[1]="Bridge Oxygen"
                elif b==3:
                    tag[0]="Oxygen attached to layer 3"
                    tag[1]="Epoxy group shifting betwen 2/3 bonds"
                else:
                    tag[0]="Between layers 2&3"

            elif l4maxz<zo<l3minz:
                if b==0:
                    tag[0]="Between layers 3&4"
                    tag[1]="Floating Oxygen atom"
                elif b==1 and bc1id_to_oxygen==1:
                    tag[0]="Between layers 3&4"
                    tag[1]="Floating C-O bond"
                elif b==1 and (c1id_to_oxygen in d["layer3_id"]):
                    tag[0]="Oxygen attached to layer 3"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==1 and (c1id_to_oxygen in d["layer4_id"]):
                    tag[0]="Oxygen attached to layer 4"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==2 and (c1id_to_oxygen in d["layer3_id"]) and (c2id_to_oxygen in d["layer3_id"]):
                    tag[0]="Oxygen attached to layer 3"
                    tag[1]="Epoxy group"
                elif b==2 and (c1id_to_oxygen in d["layer4_id"]) and (c2id_to_oxygen in d["layer4_id"]):
                    tag[0]="Oxygen attached to layer 4"
                    tag[1]="Epoxy group"
                elif b==2 and (((c1id_to_oxygen in d["layer3_id"]) and (c2id_to_oxygen in d["layer4_id"])) or ((c1id_to_oxygen in d["layer4_id"]) and (c2id_to_oxygen in d["layer3_id"]))):
                    tag[0]="Between layers 3&4"
                    tag[1]="Bridge Oxygen"
                elif b==3:
                    tag[0]="Oxygen attached to layer 4"
                    tag[1]="Epoxy group shifting betwen 2/3 bonds"
                else:
                    tag[0]="Between layers 3&4"

            elif l5maxz<zo<l4minz:
                if b==0:
                    tag[0]="Between layers 4&5"
                    tag[1]="Floating Oxygen atom"
                elif b==1 and bc1id_to_oxygen==1:
                    tag[0]="Between layers 4&5"
                    tag[1]="Floating C-O bond"
                elif b==1 and (c1id_to_oxygen in d["layer4_id"]):
                    tag[0]="Oxygen attached to layer 4"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==1 and (c1id_to_oxygen in d["layer5_id"]):
                    tag[0]="Oxygen attached to layer 5"
                    tag[1]="Oxygen attached with single bond to layer"
                elif b==2 and (c1id_to_oxygen in d["layer4_id"]) and (c2id_to_oxygen in d["layer4_id"]):
                    tag[0]="Oxygen attached to layer 4"
                    tag[1]="Epoxy group"
                elif b==2 and (c1id_to_oxygen in d["layer5_id"]) and (c2id_to_oxygen in d["layer5_id"]):
                    tag[0]="Oxygen attached to layer 5"
                    tag[1]="Epoxy group"
                elif b==2 and (((c1id_to_oxygen in d["layer4_id"]) and (c2id_to_oxygen in d["layer5_id"])) or ((c1id_to_oxygen in d["layer5_id"]) and (c2id_to_oxygen in d["layer4_id"]))):
                    tag[0]="Between layers 4&5"
                    tag[1]="Bridge Oxygen"
                elif b==3:
                    tag[0]="Oxygen attached to layer 5"
                    tag[1]="Epoxy group shifting betwen 2/3 bonds"
                else:
                    tag[0]="Between layers 4&5"

            elif zo<l5minz:
                tag[0]="Pass through all layers"
            else:
                tag[0]="Unidentified"
                print("Unidentified")

            #tag1=[1]
            #for i,(i+2),(i-38) in print_function:
            #    index1=pno.index(i)
            #    index2=pno.index(i+2)
            #    index3=pno.index(i-38)
            #    if cb1_2==3 and (bid_1[index1]==((i-41) or (i-1)) and (bid_1[index2]==((i+41) or (i+3)) and (bid_1[index3]==((i-37) or (i-39)) and (bid_2[index1]==((i-41) or (i-1)) and (bid_2[index2]==((i+41) or (i+3)) and (bid_2[index3]==((i-37) or (i-39)):
            #         tag1[0]="Crater in layer 1"   

            print(tag)
            
            f.write("%d"%eV)
            f.write(",%s"%tag[0])
            f.write(",%s"%tag[1])
            f.write(",%s"%tag[2])
            #f.write(",%s"%tag[3])
            #f.write(",%.3f"%zo)
            f.write(",%.3f"%heightofO)
            #f.write(",%.3f"%Xposo)
            #f.write(",%.3f"%Yposo)
            f.write("\n")
            #print(tag1)
f.close()


