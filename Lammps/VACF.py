import matplotlib.pyplot as plt
from fractions import Fraction
import numpy as np
import math

# 375 for 3000 & 5000 run
dt=0.25
Dumpstepsize=25
Totalrun = 30000
Diffcoeffwindow=[]
for u in range(5,6):
	Window = (u+1)*2
	Startrun=1125
	#Diffcoeff=[]


	for t in range((Startrun/Dumpstepsize)+1,((Totalrun+Dumpstepsize)/Dumpstepsize)):
		run_i=t*Dumpstepsize
		run = str(run_i)
		item = "TIMESTEP"

		f=open("diffusion925K_1.txt","r")
		file=f.readlines()

		#OBTAIN VELOCITIES OF O ATOM IN FILE
		vx=[]
		vy=[]
		vz=[]
		count=0
		for line in file:
			count += 1
		linesno=count

		check1=[]
		check2=[]

		for lines in file:
			check1.append(lines.split()[0])
			try:
				check2.append(lines.split()[1])
			except IndexError:
				check2.append(" ")
			continue

		for i in range(0,(linesno-1)):
			if check1[i+1] == str(Startrun) and check2[i]==item:
				h=i+2

		for i in range(0,(linesno-1)):
			if check1[i+1] == run and check2[i]==item:
				p=i+2

		start=h+8
		end=p+8
		#print(start)
		x=file[(start-1):(end):10]
		for line in x:
		    vx.append(float(line.split()[1]))
		    vy.append(float(line.split()[2]))
		    vz.append(float(line.split()[3]))

		f.close()

	#SET SIZE OF ARRAYS
	#vacfb=np.array([0]*(len(vx)-Window),dtype='f')
	indices =list(range(1,(len(vx))))
	indices1 =list(range(0,(len(vx))-1))
	indices2 =list(range(0,(len(vx))-2))
	indices3 =list(range(1,(len(vx))-2))
	#Tau = [l*Dumpstepsize*dt for l in indices]
	#print(tstep)
	#print(len(tstep))

	#print(len(x))
	#print(Window)
	#print(int(floor(Fraction(int(len(x)), int(Window)))))
	#torigin=list(range(0,Window*(int(floor(Fraction(int(len(x)), int(Window))))+1),Window))
	#torigin = [0]*(len(indices))
	#print(torig[0])

	deltat =list(range(0,(len(vx)-1),Window))
	addx = [0]*len(vx)

	vx=vx+addx
	vy=vy+addx
	vz=vz+addx

	#torigin = [[p for p in indices if (p+l)<(len(indices)+1)] for l in indices]
	#print(torigin[1153])
	#print(torig[1150])
	torigin=[0]*len(indices1)
	for l in indices1:
			torigin[l] = [ele for ele in deltat if ele<=l]

	#print(torig)
	#originsfort = [[p if p<L for l in tstep] for l in tstep]
	#divide = [int(len(torig[l])) for l in torig]
	#print(x[1155])
	sumx = [np.sum([(vx[l+p]*vx[l]) for p in torigin[l]]) for l in indices1]
	sumy = [np.sum([(vy[l+p]*vy[l]) for p in torigin[l]]) for l in indices1]
	sumz = [np.sum([(vz[l+p]*vz[l]) for p in torigin[l]]) for l in indices1]

	#sumx = [np.sum([(vx[l+p]*vx[l]) for p in torigin[l]]) for l in indices1]
	#sumy = [np.sum([(vx[p]*vy[l+p]) for p in torigin[l]]) for l in indices1]
	#sumz = [np.sum([(vz[p]*vz[l+p]) for p in torigin[l]]) for l in indices1]

	#print(sumx)
	#print(tstep1)

	vacfx=[float(1/ float(len(torigin[l])))*float(sumx[l]) for l in indices2]
	vacfy=[float(1/ float(len(torigin[l])))*float(sumy[l]) for l in indices2]
	vacfz=[float(1/ float(len(torigin[l])))*float(sumz[l]) for l in indices2]

	Time2=[l*Dumpstepsize*dt for l in indices]

	VACF = [vacfx[l]+vacfy[l]+vacfz[l] for l in indices2]

	#Diffcoeff= 
	#print(len(x))
	#print(Window)
	#print(int(floor(Fraction(int(len(x)), int(Window)))))
	#torigin=list(range(0,Window*(int(floor(Fraction(int(len(x)), int(Window))))+1),Window))
	#timedoubleprime = [0]*(len(indices))
	#print(len(timedoubleprime))
	#print(torig[0])

	#timedoubleprime = [[p for p in indices1 if p<(l-1)] for l in indices1]
	#print(timedoubleprime[30])
	#print(timedoubleprime[0])
	#print(timedoubleprime[1114])

	#sumx = [np.sum([np.sum([(vx[l-k]*vx[p]) for k in timedoubleprime[l]]) for p in torig[l]]) for l in indices1]
	#vacfx=[float(float(1)/ float(len(torig[l])))*float(sumx[l]) for l in indices1]
	#print(vacfx)
	Time=[l*Dumpstepsize*dt*0.001 for l in indices1]
	Diffcoeff=[float(Fraction(int(1), int(3)))*np.sum(VACF[0:l])*Dumpstepsize*dt*np.exp(-5) for l in indices2]
	#vacf=[l*np.exp(10) for l in vacffs]

	print(Diffcoeff[len(Diffcoeff)-1])
	#print(len(vacffs))
	Diffcoeffwindow.append([l for l in Diffcoeff])
	#PLOT GRAPHS
	print(Diffcoeffwindow[0])


	plt.subplot(212)
	plt.plot(Time[0:len(Diffcoeff)], Diffcoeffwindow[u], label="Window=%d"%Window)
	plt.title('Diffusion coefficient vs. Time at 700K')
	plt.xlabel('Picoseconds [ps]')
	plt.ylabel('Diffusion Coefficient [m^2/s]')
	plt.annotate('%0.6f' % Diffcoeff[len(Diffcoeff)-1].max(), xy=(1, Diffcoeff[len(Diffcoeff)-1].max()), xytext=(8, 0), 
	                 xycoords=('axes fraction', 'data'), textcoords='offset points')
	plt.legend(loc='centerleft',bbox_to_anchor=(1,0.5))
plt.subplot(211)
plt.plot(Time2[0:len(VACF)], VACF, label="%d"%Startrun)
plt.title('VACF vs. Time at 700K')	
plt.xlabel('femtoseconds [fmts]')
plt.ylabel('VACF [(A/fmts)^2]')	
plt.show()
	#plt.legend(loc='centerleft',bbox_to_anchor=(1,0.5))
