import matplotlib.pyplot as plt
import numpy as np
def plotcmMSD():
	with open("cmMSD","r") as c:
		with ("XDATCAR","r") as x:
			L=np.array([[0 for i in range(3)]for j in range(3)])
			for i in range(7):
				line=f.readline()
				if(i==2):
					L[0]=list(map(float,line.split()))
				if(i==3):
					L[1]=list(map(float,line.split()))
				if(i==4):
					L[2]=list(map(float,line.split()))
		temp=np.empty((0,4), float)
		for line in c:
			temp=np.append(temp,np.array([list(map(float,line.split()))]),axis=0)
		tc=temp.T
		N=len(temp)
		plt.plot(tc[0],tc[1],label="x")
		plt.plot(tc[0],tc[2],label="y")
		plt.plot(tc[0],tc[3],label="z")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$r(Å))$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.legend()
		plt.savefig('cmMSDxyz.eps', dpi=150)
		plt.savefig('cmMSDxyz.png', dpi=150)
		plt.tight_layout()
		plt.show()
		
		print("Input Lattice param")
		temp=tc[1]*L[0]+tc[2]*L[1]+tc[3]*L[2]
		
		MSD=np.sum([[tc[j][i]*tc[j][i] for i in range(N)]for j in range(1,4)],axis=0)
		MSD=np.sum([[temp[j][i]*temp[j][i] for i in range(N)]for j in range(1,4)],axis=0)
		#aboveMSDis error
		MSD=MSD*L**2
		plt.plot(tc[0],MSD,label="cmMSD")
		plt.legend(loc="lower right")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('cmMSDliner.eps', dpi=150)
		plt.savefig('cmMSDliner.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
		plt.plot(tc[0],MSD,label="cmMSD")
		plt.legend(loc="lower right")
		plt.yscale("log")
		plt.xscale("log")
		plt.grid(which="both")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('cmMSDlog.eps', dpi=150)
		plt.savefig('cmMSDlog.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
def plotMSD():
	with open("cmSiMSD","r") as cS,open("cmOMSD","r") as cO:
		temp=np.empty((0,2), float)
		for line in cS:
			temp=np.append(temp,np.array([list(map(float,line.split()))]),axis=0)
		ts=temp.T
		temp=np.empty((0,2), float)
		for line in cO:
			temp=np.append(temp,np.array([list(map(float,line.split()))]),axis=0)
		to=temp.T
		plt.plot(ts[0],ts[1],label="Si")
		plt.plot(to[0],to[1],label="O")
		plt.legend(loc='lower right')
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('SiOMSDliner.eps', dpi=150)
		plt.savefig('SiOMSDliner.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		plt.plot(ts[0],ts[1],label="Si")
		plt.plot(to[0],to[1],label="O")
		plt.legend(loc='lower right')
		plt.yscale("log")
		plt.xscale("log")
		plt.grid(which="both")
		plt.xlabel("Time(fs)", fontsize=16)
		plt.ylabel("$MSD(Å^{2})$", fontsize=16)
		plt.tick_params(labelsize=16)
		plt.savefig('SiOMSDlog.eps', dpi=150)
		plt.savefig('SiOMSDlog.png', dpi=150)
		plt.tight_layout()
		plt.clf()
		
if __name__ =="__main__":
	#plotcmMSD()
	plotMSD()