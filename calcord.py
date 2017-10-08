import numpy as np
import matplotlib.pyplot as plt
print("This program only works on a-quartz+SiO")
def calcoord():
	with open("XDATCAR","r") as f,open("cord.dat","w") as g:
		
		f.readline()
		f.readline()
		L=float(f.readline().split()[0])
		f.readline()
		f.readline()
		atoms=(f.readline().split())
		atomnum=list(map(int,f.readline().split()))
		atomlist=[]
		for (a,b) in zip(atoms,atomnum):
			atomlist.append([a,b])
		print(atomnum)
		sumnum=sum(atomnum)
		P=np.empty((0,3), float)
		Si=np.empty((0,3), float)
		O=np.empty((0,3), float)
		N=np.empty((0,3), float)
		coordnum=[]
		linecount=0
		time=0
		short=2.4
		long=3
		flag=0
		for line in f:
			if(linecount==0 or linecount%(sumnum+1)==0):
				if(flag==1):
					coordnum.append(calculate_coord_count(P,O,time,short,long,L))	#calprocess
					g.write(str(coordnum[0][0])+" "+str(coordnum[0][1])+" "+str(coordnum[0][2])+"\n")
				
				P=np.empty((0,3), float) #P=addSi
				Si=np.empty((0,3), float)
				O=np.empty((0,3), float)
				N=np.empty((0,3), float)
				linecount+=1
				
				continue
			elif(linecount%(sumnum+1)<atomlist[0][1]+1):
				P=np.append(P,np.array([list(map(float,line.split()))]),axis=0)
				linecount+=1
				continue
			elif(linecount%(sumnum+1)<atomlist[0][1]+atomlist[1][1]+1):
				Si=np.append(Si,np.array([list(map(float,line.split()))]),axis=0)
				linecount+=1
				continue
			elif(linecount%(sumnum+1)<atomlist[0][1]+atomlist[1][1]+atomlist[2][1]+1):
				O=np.append(O,np.array([list(map(float,line.split()))]),axis=0)
				linecount+=1
				continue
			else:
				N=np.append(N,np.array([list(map(float,line.split()))]),axis=0)
				linecount+=1
				time+=0.001
				flag=1
				continue
				
		plt.plot(coordnum[0][0],coordnum[0][1],label="short")
		plt.plot(coordnum[0][0],coordnum[0][2],label="long")
		plt.legend()
		plt.show()
def calculate_coord_count(atom1,atom2,T,s,l,a):
	#This calculate short and long distance coordination number of between atom1 and atom2
	#roundを使って周期的境界条件を考慮する
	shortnum,longnum=0,0
	for el in atom2:
		sub=np.array([round(el[0]-atom1[0][0]),round(el[1]-atom1[0][1]),round(el[2]-atom1[0][2])])
		dis=(np.dot(sub,sub))**(1/2)*a
		if dis<=l:
			longnum+=1
			if dis<=s:
				shortnum+=1
			
	return T,shortnum,longnum
calcoord()