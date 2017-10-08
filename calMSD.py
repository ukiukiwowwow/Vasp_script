import copy
import re
import numpy as np
import sys
def MSD(flag=True):
	with open("XDATCAR","r") as f,open("cmMSD","w") as cM,open("MSD","w") as M,open("r.dat","w") as R:
		
		L=np.array([np.zeros(3)]*3)
		for i in range(7):
			line=f.readline()
			if(i==1):
				mult=float(line)
			if(i==2):
				L[0]=list(map(float,line.split()))
			if(i==3):
				L[1]=list(map(float,line.split()))
			if(i==4):
				L[2]=list(map(float,line.split()))
			if(i==5):
				atoms=line.split()
			if(i==6):
				atomnum=list(map(int,line.split()))
		L=L*mult
		print(L)
		#方針を変えて一つのMSDファイルにまとめる Si O/time SiMSD OMSDのようにする？ 
		if flag:
			EA=open("eachMSD","w")
			for i in range(len(atoms)):
				for j in range(atomnum[i]):
					EA.write(str(atoms[i])+","+str(j+1)+" ")
			EA.write("\n")
		
		for i in atoms:
			M.write(str(i)+" ")
			
		#exec('cm{}MSD=open("cm"+str(i),"w")'.format(i))#cmSiMSDの後cmOMSDを作ると前のcmSiMSDのような変数が破棄される？
		M.write("\n")
		sumatom=sum(atomnum)
		cur=np.array([np.zeros(3)]*sumatom)
		pre=np.array([np.zeros(3)]*sumatom)
		r=np.array([np.zeros(3)]*sumatom)
		cmMSD=np.zeros(3)
		d=np.zeros(3)
		#Totalstep=int(input())
		time=-0.001
		#Tlines=(sumatom+1)*(Totalstep)
		linecount=0
		echo=0
		#msditerに渡す？
		while True:
			line=f.readline()
			if(line==""):
				print("OK?")
				break
			if(linecount%(sumatom+1)==0):
				#MSD計算
				#msd,cmmsd=calmsd(cur,pre,atom,atomnum)
				if(pre.any()==True):
					for an in range(sumatom):
						for dim in range(3):
							d[dim]=cur[an][dim]-pre[an][dim]
							d[dim]-=round(d[dim])
							r[an][dim]+=d[dim]
					cmMSD=np.sum(r,axis=0)/sumatom
					cM.write(str(time)+" "+str(cmMSD[0])+" "+str(cmMSD[1])+" "+str(cmMSD[2])+"\n")
					#Sir=(Six-cmMSD[0])*L[0]+(Siy-cmMSD[1])*L[1]+(Siz-cmMSD[2])*L[2]
					#Simsd=1/(Sinum)*np.dot(Sir,Sir)
					#ps.write(str(time)+" "+str(Simsd)+"\n")
					atomr=(r-cmMSD)[:,0].reshape(1,len(r)).T*L[0]+(r-cmMSD)[:,1].reshape(1,len(r)).T*L[1]+(r-cmMSD)[:,2].reshape(1,len(r)).T*L[2]#3*N次元の物を3*N次元にする
					#atomrがおかしい r,cmMSDは間違いなく同一である。Lがおかしかったです
					temp=0;M.write(str(time)+" ")
					if flag:
						EA.write(str(time)+" ")
						for i in range(sumatom):
							eachmsd=(np.linalg.norm(atomr[i])**2)
							EA.write(str(eachmsd)+" ")
						EA.write("\n")
					for i in range(len(atoms)):
						Latom=atomr[temp:temp+atomnum[i]]#-1?
						"""
						if(i==0):
							R.write(str(Latom)+"\n")
						"""
						msd=(np.linalg.norm(Latom)**2)/atomnum[i]
						M.write(str(msd)+" ")#exec('cm{}MSD.write(str(time)+" "+str(msd)+"\n")'.format(i))
						temp=atomnum[i]#-1?
					M.write("\n")
				linecount=1
				time+=0.001
				echo+=1
				if(echo%100==0):
					print(echo)
				pre=copy.deepcopy(cur)
				cur=np.array([np.zeros(3)]*sumatom)
				continue
			else:
				
				cur[linecount-1]=list(map(float,line.split()))
				linecount+=1
				continue
		"""	
		for i in atoms:
			exec('cm{}.close()'.format(i))
		"""
		if flag:
			EA.close()
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
	MSD(flag)