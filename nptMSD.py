import copy
import re
import numpy as np
import sys
import IPython
def nptMSD(flag=True):
	"""
	エラー原因
	読み込み範囲を間違えているため msdの値が変なことになっているのでは
	cur[??] ??がおかしいですね
	"""
	with open("XDATCAR","r") as f,open("cmMSD","w") as cM,open("MSD","w") as M,open("r.dat","w") as R:
		
		L=np.array([np.zeros(3)]*3)
		li=8
		isFirst=True
		time=0.
		d=np.zeros(3);count=0
		cmMSD=np.zeros(3)
		for index,line in enumerate(f):
			if(isFirst):
				if(index==0):
					continue
				elif(index==1):
					mult=float(line)
					continue
				elif(index==2):
					L[0]=list(map(float,line.split()))
					continue
				elif(index==3):
					L[1]=list(map(float,line.split()))
					continue
				elif(index==4):
					L[2]=list(map(float,line.split()))
					continue
				elif(index==5):
					atomname=line.split()
					continue
				elif(index==6):
					atomnum=list(map(int,line.split()))
					sumatom=sum(atomnum)
					#isFirst=False
					cur=np.array([np.zeros(3)]*sumatom)
					pre=np.array([np.zeros(3)]*sumatom)
					r=np.array([np.zeros(3)]*sumatom)
					continue
				elif(index==7):
					continue
				else:
					cur[index-8]=list(map(float,line.split()));print("First {0}".format(index-8))
					#cur=np.append(cur,list(map(float,line.split())))
					if(index==(sumatom+7) and len(pre)==sumatom):
						isFirst=False;print("?")
			else:
				if(index%(li+sumatom)==0):
					if(pre.any()==True):
						r+=cur-pre
						cmMSD=np.sum(r,axis=0)/sumatom
						cM.write("{0} {1} {2} {3}\n".format(time,cmMSD[0],cmMSD[1],cmMSD[2]))
						atomr=(r-cmMSD)[:,0].reshape(1,len(r)).T*L[0]+(r-cmMSD)[:,1].reshape(1,len(r)).T*L[1]+(r-cmMSD)[:,2].reshape(1,len(r)).T*L[2]
						M.write("{0} ".format(time))
						temp=0
						for i in range(len(atomname)):
							Latom=atomr[temp:temp+atomnum[i]]#-1?
							msd=(np.linalg.norm(Latom)**2)/atomnum[i]
							M.write(str(msd)+" ")
							temp=atomnum[i]
						M.write("\n")
						
					pre=copy.deepcopy(cur)
					cur=np.array([np.zeros(3)]*sumatom)
					time+=0.001 # 1fs unfixed
					count=0
					continue
				elif(index%(li+sumatom)==2):
					k=line.split()
					L[0]=list(map(float,line.split()))
					continue
				elif(index%(li+sumatom)==3):
					L[1]=list(map(float,line.split()))
					continue
				elif(index%(li+sumatom)==4):
					L[2]=list(map(float,line.split()))
					continue
				elif(index%(li+sumatom)==5 or index%(li+sumatom)==6 or index%(li+sumatom)==7 or index%(li+sumatom)==1):
					continue
				else:
					cur[count]=list(map(float,line.split()));print("notFirst {0}".format(count))
					count+=1
					continue

		r+=cur-pre
		cmMSD=np.sum(r,axis=0)/sumatom
		cM.write("{0} {1} {2} {3}\n".format(time,cmMSD[0],cmMSD[1],cmMSD[2]))
		atomr=((r-cmMSD)[:,0].reshape(1,len(r)).T*L[0]+(r-cmMSD)[:,1].reshape(1,len(r)).T*L[1]+(r-cmMSD)[:,2].reshape(1,len(r)).T*L[2])
		M.write("{0} ".format(time))
		temp=0
		#IPython.embed()
		for i in range(len(atomname)):
			Latom=atomr[temp:temp+atomnum[i]]#-1?
			msd=(np.linalg.norm(Latom)**2)/atomnum[i]
			M.write(str(msd)+" ")
			temp=atomnum[i]
			M.write("\n")
		IPython.embed()
if __name__ =="__main__":
	print("python calMSD.py -eでeachMSD")
	flag=True
	try:
		t=sys.argv[1]
		if(t=="-e"):
			flag=True
	except:
		print("If you want to calculate each atom,strike \"y\" key")
		if input("[y] or anything else >> ")[0].lower=="y":
			flag=True
		else:
			flag=False
	nptMSD(flag)
