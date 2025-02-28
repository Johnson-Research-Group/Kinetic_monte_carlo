import matplotlib.pyplot as plt
import numpy as np

# 375 for 3000 & 5000 run
dt=0.25
Dumpstepsize=25
Totalrun = 60100
nsamp = 200
dtime=dt*nsamp
for u in range(45,46):
	Startrun=u*25 #1125
	Diffcoeff=[]
	Time=[]


	for t in range((Startrun/Dumpstepsize)+1,((Totalrun+Dumpstepsize)/Dumpstepsize)):
		run_i=t*Dumpstepsize
		run = str(run_i)
		item = "TIMESTEP"
	    
		ti=t-((Startrun/Dumpstepsize)+1)

		f=open("diffusion300K_2.txt","r")
		file=f.readlines()

		#LOCATING VELOCITIES OF O ATOM IN FILE
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

		print(run_i)
		print(h)
		print(p)
		start=h+8
		end=p+8
		print(start)
		print(end)
		x=file[(start-1):(end):10]
		for line in x:
		    vx.append(float(line.split()[1]))
		    vy.append(float(line.split()[2]))
		    vz.append(float(line.split()[3]))

		f.close()

		#Initial Velocities
		vx0 = vx[0]
		vy0 = vy[0]
		vz0 = vz[0]

		#VELOCITIY AUTO-CORRELATION
		print(vx)
		#Make empty array of the length vx[]
		vacf=[None]*(len(vx))
		
		for i in range(0,(len(vx))):
			vacf[i] = vx[i]*vx0 + vy[i]*vy0 + vz[i]*vz0

		a = vacf[0]
		b = vacf[(len(vx)-1)]

		vacf[0] = 0.5*a
		vacf[(len(vx)-1)] = 0.5*b

		trap_sum = np.sum(vacf)
	    
	    #DIFFUSION COEFFICIENT
		Diffcoeff.append((dt*trap_sum)*np.exp(-5))

		print(run_i)
		print("Diffusion coefficient = %.6f m^2/s" %(Diffcoeff[ti]))
	    
		Time.append((run_i-Startrun)*0.25*np.exp(-15))
		

	print("HERE IT IS: %d" %u)

	#PLOT GRAPHS
	plt.plot(Time, Diffcoeff, label="%d"%Startrun)
plt.title('Diffusion coefficient vs. Runtime')
#plt.legend(loc='centerleft',bbox_to_anchor=(1,0.5))
#plt.plot(Runtime, Diffcoeff)
plt.xlabel('Seconds [s]')
plt.ylabel('Diffusion Coefficient [m^2/s]')
plt.show()
#print(vz)

