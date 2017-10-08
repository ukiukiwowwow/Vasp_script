import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
#もしnewMSDファイルを使うなら
def newdiff(flag=1):
	with open("MSD","r") as f,open("subdiff.dat","w") as sd: 
		atom=f.readline().split()
		atomname="".join(atom)
		sumDiff=np.array([.0]*len(atom))
		for i in atom:
			sd.write(str(i)+" ")
		sd.write("\n")
		if(flag==1):
			temp=np.array([])
			temp2=np.array([])
		linecount=0
		count=0
		for line in f:
			linecount+=1
			l=list(map(float,line.split()))
			time,atoms=l[0],np.array(l[1:])
			temp2=np.append(temp2,np.array([list(map(float,line.split()))]))
			if(flag==1):
				temp=np.append(temp,time)
				temp=np.append(temp,np.array(atoms)/(6*time)*10)#あとでreshapeをする
			if(linecount%5000==0):
				print("{0},{1}".format(time,atoms/(6*time)*10))
				sd.write(str(time)+" ")
				for i in range(len(atom)):
					sd.write(str(atoms[i]/(6*time)*10)+" ")
				sd.write("\n")
			if(linecount>10000):
				sumDiff+=atoms/(6*time)*10
				count+=1
		print(time,atoms/(6*time)*10)
		print("Mean D are {0}".format(sumDiff/count))
		sd.write(str(time)+" ")
		for i in range(len(atom)):
			sd.write(str(atoms[i]/(6*time)*10)+" ")
		sd.write("\n"+"Mean D are "+str(sumDiff/count)+"\n")
		if(flag==1):
			temp=temp.reshape(len(temp)/(1+len(atom)),(1+len(atom)))
			temp2=temp2.reshape(len(temp2)/(1+len(atom)),(1+len(atom)))
			D=temp.T
			MSD=temp2.T
		print(len(D),len(D[0]))
		if(flag==1):
			for i in range(1,len(atoms)+1):
				plt.plot(D[0],D[i],label=str(atom[i-1]))
			plt.legend(loc=1)
			plt.xlabel("Time(fs)", fontsize=16)
			plt.ylabel("$D(10^{-9}m^{2}/s)$", fontsize=16)
			plt.tick_params(labelsize=16)
			plt.savefig('Dliner.eps', dpi=150)
			plt.savefig('Dliner.png', dpi=150)
			plt.tight_layout()
			plt.clf()
			
			for i in range(1,len(atoms)+1):
				plt.plot(D[0],D[i],label=str(atom[i-1]))
			plt.legend(loc=1)
			plt.yscale("log")
			plt.xscale("log")
			plt.grid(which="both")
			plt.xlabel("Time(ps)", fontsize=16)
			plt.ylabel("$D(10^{-9}m^{2}/s)$", fontsize=16)
			plt.tick_params(labelsize=16)
			plt.savefig('Dlog.eps', dpi=150)
			plt.savefig('Dlog.png', dpi=150)
			plt.tight_layout()
			plt.clf()
			
			for i in range(1,len(atoms)+1):
				plt.plot(MSD[0],MSD[i],label=str(atom[i-1]))
			plt.legend(loc='lower right')
			plt.xlabel("Time(ps)", fontsize=16)
			plt.ylabel("$MSD(Å^{2})$", fontsize=16)
			plt.tick_params(labelsize=16)
			
			plt.savefig('MSDliner.eps', dpi=150)
			plt.savefig('MSDliner.png', dpi=150)
			plt.tight_layout()
			plt.clf()
		
			for i in range(1,len(atoms)+1):
				plt.plot(MSD[0],MSD[i],label=str(atom[i-1]))
			plt.legend(loc='lower right')
			plt.yscale("log")
			plt.xscale("log")
			plt.grid(which="both")
			plt.xlabel("Time(ps)", fontsize=16)
			plt.ylabel("$MSD(Å^{2})$", fontsize=16)
			plt.tick_params(labelsize=16)
			plt.savefig('MSDlog.eps', dpi=150)
			plt.savefig('MSDlog.png', dpi=150)
			plt.tight_layout()
			plt.clf()
			
def cmdiff(flag=1):
	with open("cmMSD","r") as cM,open("cmdiff.dat","w") as cd:
		print("cm-selfdiffusion")
		cd.write("cm-selfdiffusion\n")
		with open("XDATCAR","r") as x:
			L=np.array([[0.1 for i in range(3)]for j in range(3)])
			for i in range(7):
				line=x.readline()
				if(i==1):
					mult=float(line)
				if(i==2):
					L[0]=list(map(float,line.split()))
				if(i==3):
					L[1]=list(map(float,line.split()))
				if(i==4):
					L[2]=list(map(float,line.split()))
			L=mult*L
		temp=np.empty((0,2), float)
		T=np.array([])
		X=np.array([])
		Y=np.array([])
		Z=np.array([])
		MSD=np.array([])
		linecount=0
		for line in cM:
			linecount+=1
			time,x,y,z=list(map(float,line.split()))
			T=np.append(T,time)
			X=np.append(X,x)
			Y=np.append(Y,y)
			Z=np.append(Z,z)
			r=x*L[0]+y*L[1]+z*L[2]
			r=np.dot(r,r)
			MSD=np.append(MSD,r)
			#r=np.dot(x,x)+np.dot(y,y)+np.dot(z,z)
			#r=r*L**2
			temp=np.append(temp,np.array([[time,(r/(6*time)*10)]]),axis=0)
			if(linecount%5000==0):
				print("{0},{1}".format(time,r/(6*time)*10))
				cd.write(str(time)+" "+str(r/(6*time)*10)+"\n")
		print(time,r/(6*time)*10)
		
		
		plt.plot(T,X,label="x")
		plt.plot(T,Y,label="y")
		plt.plot(T,Z,label="z")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$r(Å))$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.legend()
		plt.savefig('cmMSDxyz.eps', dpi=150)
		plt.savefig('cmMSDxyz.png', dpi=150)
		plt.tight_layout()
		plt.show()
		plt.clf()
		
		plt.plot(T,MSD,label="cmMSD")
		plt.legend(loc="lower right")
		plt.xlabel("Time(ps)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('cmMSDliner.eps', dpi=150)
		plt.savefig('cmMSDliner.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
		plt.plot(T,MSD,label="cmMSD")
		plt.legend(loc="lower right")
		plt.yscale("log")
		plt.xscale("log")
		plt.grid(which="both")
		plt.xlabel("Time(ps)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('cmMSDlog.eps', dpi=150)
		plt.savefig('cmMSDlog.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
		
		Dc=temp.T
		plt.plot(Dc[0],Dc[1],label="Dcm")
		plt.legend(loc=1)
		plt.xlabel("Time(ps)", fontsize=16)
		plt.ylabel("$D(10^{-9}m^{2}/s)$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('Dcmliner.eps', dpi=150)
		plt.savefig('Dcmliner.png', dpi=150)
		plt.tight_layout()
		plt.clf()
			
		plt.plot(Dc[0],Dc[1],label="Dcm")
		plt.legend(loc=1)
		plt.yscale("log")
		plt.xscale("log")
		plt.grid(which="both")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$D(10^{-9}m^{2}/s)$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('Dcmlog.eps', dpi=150)
		plt.savefig('Dcmlog.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
		
def eachdiff():
	with open("../eachMSD","r") as f:
		
		atom=f.readline().split()
		atomname=[]
		for item in atom:
			atomname.append([item[0],item[1]])
		temp=np.array([])
		linecount=0
		for line in f:
			linecount+=1
			l=list(map(float,line.split()))
			time,atoms=l[0],np.array(l[1:])
			temp=np.append(temp,np.array(l))
			if(linecount%5000==0):
				print(linecount)
		
		temp=temp.reshape(len(temp)/(1+len(atom)),(1+len(atom))) #(行,列)
		MSD=temp.T
		for i in range(1,len(atoms)+1):
			plt.plot(MSD[0],MSD[i],label=str(atom[i-1]))
			plt.legend(loc='lower right')
			plt.xlabel("Time(ps)", fontsize=16)
			plt.ylabel("$MSD(Å^{2})$", fontsize=16)
			plt.tick_params(labelsize=16)
			plt.show()
			plt.savefig(str(atom[i-1].split(",")[0])+"MSD"+str(i)+"liner.eps", dpi=150)
			plt.savefig(str(atom[i-1].split(",")[0])+"liner.png", dpi=150)
			plt.tight_layout()
			plt.clf()
		
		for i in range(1,len(atoms)+1):
			plt.plot(MSD[0],MSD[i],label=str(atom[i-1]))
			plt.legend(loc='lower right')
			plt.yscale("log")
			plt.xscale("log")
			plt.grid(which="both")
			plt.xlabel("Time(ps)", fontsize=16)
			plt.ylabel("$MSD(Å^{2})$", fontsize=16)
			plt.tick_params(labelsize=16)
			plt.savefig(str(atom[i-1].split(",")[0])+"MSD"+str(i)+"log.eps", dpi=150)
			plt.savefig(str(atom[i-1].split(",")[0])+"MSD"+str(i)+"log.png", dpi=150)
			plt.tight_layout()
			plt.clf()
		
if __name__=="__main__":
	#print("Do you want to plot D,type(y)")
	#if re.compile("y",re.IGNORECASE).match(input().split()[0]) != None:
	flag=1
	#else:
	#	flag=0
	newdiff(flag)
	cmdiff(flag)
	if(os.path.exists("./eachMSD")):
		try:
			os.mkdir("./eachatomplot");os.chdir("./eachatomplot")
		except:
			os.chdir("./eachatomplot")
		print(os.getcwd())
			
		eachdiff()
	