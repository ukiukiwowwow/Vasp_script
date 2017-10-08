import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import time
import copy
def cord_time_change(sisi,sio,dx):
	"""sisiはSi-Siの結合距離 sioはSi-Oの結合距離
	方針として各Si[i]ごとに結合パターンSi:O=[x,y]を求める
	それを時間との表にして順次出力する。(原子ごとに異なるファイルにする？)
	プロットは別にやる(理由空間計算量大)
	入れ替わったかどうかを判断するのは難しく入れ替わりだけを対象としているのでそれ以外は捕らえられない
	
	
	Si原子ごとに周囲のOの数とSiの数を記録する
	
	
	
	各原子の変位量を見る
	ある原子をi番目としてある時刻をt0として
	t0でのiの位置をxi0 t0+10fsの位置をxi10でこれが一定以上大きくなると
	つまりif(xi10-xi0)のときのi原子の配位数の変化を見る
	t0での配位数ci0 t+10fsでの配位数ci10 このときのci0とci10を記録する。
	そしてxi10 xi0 ci10 ci0を更新して次のステップへ
	"""
	timecount=0
	flag=False
	print(os.getcwd())
	shutil.copyfile("../XDATCAR","./XDATCAR")
	print("copy")
	with open("XDATCAR","r") as f,open("Sitimecord.dat","w") as Sc,open("precord.dat","w") as pc,open("curcord.dat","w") as cc:
		L=np.array([np.zeros(3)]*3)
		for i in range(5):
			line=f.readline()
			if(i==2):
				L[0]=list(map(float,line.split()))
				continue
			if(i==3):
				L[1]=list(map(float,line.split()))
				continue
			if(i==4):
				L[2]=list(map(float,line.split()))
				continue
		print(L)
		atomname=np.array(f.readline().split())
		atomnum=np.array(list(map(int,f.readline().split())))
		sumatom=np.sum(atomnum)
		#x0=np.array(np.zeros(3*atomnum[0])).reshape(atomnum,3);x10=np.array(np.zeros(3*atomnum[0])).reshape(atomnum,3)
		#csi0=np.array(np.zeros(atomnum[0]));csi10=np.array(np.zeros(atomnum[0]))
		#co0=np.array(np.zeros(atomnum[0]));co10=np.array(np.zeros(atomnum[0]))
		cSiSit=np.array([])
		cSiOt=np.array([])
		xSit=np.array([])
		"""
		x10-x0>r0 1.6Å程度
		t10-t0>10fs
		
		"""
		atom=np.array([np.zeros(3)]*sumatom)
		linecount=0
		for index in range(atomnum[0]):
			Sc.write("Si"+str(index+1)+" ")
			cc.write("Si"+str(index+1)+" ")
			pc.write("Si"+str(index+1)+" ")
		cc.write("\n")
		pc.write("\n")
		Sc.write("\n")
		while True:
			line=f.readline()
			if(line==""):
				print("OK?")
				break
			if(linecount%(sumatom+1)==0):
				#Siの各結合数を数える
				#Si-O
				SiOnum=np.zeros(atomnum[0])
				SiSinum=np.zeros(atomnum[0])
				if atom.any():
					for Siindex,Si in enumerate(atom[:atomnum[0]]):
						Ocount=0
						for O in atom[atomnum[0]:]:
							fracdis=Si-O
							#0.5以上の距離を持つ場合周期的境界条件の補正を加える
							for dim in range(3):
								fracdis[dim]-=round(fracdis[dim])
							realdis=L[0]*fracdis[0]+L[1]*fracdis[1]+L[2]*fracdis[2]
							#Lをかけよう
							if(np.linalg.norm(realdis)<=sio):
								Ocount+=1
						SiOnum[Siindex]=int(Ocount)
					for Siindex,Si1 in enumerate(atom[:atomnum[0]]):
						Sicount=0
						for Si2 in atom[:atomnum[0]]:
							fracdis=Si1-Si2
							#0.5以上の距離を持つばあい周期的境界条件の補正を加える
							for dim in range(3):
								fracdis[dim]-=round(fracdis[dim])
							realdis=L[0]*fracdis[0]+L[1]*fracdis[1]+L[2]*fracdis[2]
							#Lをかけよう
							if(np.linalg.norm(realdis)<=sisi):
								Sicount+=1
						SiSinum[Siindex]=int(Sicount-1)
						
					for index in range(atomnum[0]):
						Sc.write(str(SiOnum[index])+","+str(SiSinum[index])+" ")
						if(timecount>20):
							t0=xSit[:atomnum[0]*3].reshape(atomnum[0],3)
							#print(t0)
							dis=Si-t0[index]
							for dim in range(3):
								dis[dim]-=round(dis[dim])
							realdis=L[0]*dis[0]+L[1]*dis[1]+L[2]*dis[2]
							if(np.linalg.norm(realdis)>dx):#移動したと考えるsioパラメータとは別に与える可能性あり
								#print(index,timecount)
								cc.write(str(index+1)+","+str(SiOnum[index])+","+str(SiSinum[index])+" ")
								tSiO0=cSiOt[:atomnum[0]]
								tSiSi0=cSiSit[:atomnum[0]]
								pc.write(str(index+1)+","+str(tSiO0[index])+","+str(tSiSi0[index])+" ")
								flag=True
					Sc.write("\n")
					if flag:
						cc.write("\n")
						pc.write("\n")
						flag=False
					if(timecount<=10):
						cSiOt=np.append(cSiOt,SiOnum) #[t0,t1,t2,t3,t4...t10]
						cSiSit=np.append(cSiOt,SiSinum)
						xSit=np.append(xSit,atom[:atomnum[0]]) #間違っているよ 一番新しいt10が[0]にはいるよ 
					else:#治った？
						cSiOt=np.append(cSiOt[atomnum[0]:],SiOnum)#[t11,t1,t2]->[t1,..t11]
						cSiSit=np.append(cSiSit[atomnum[0]:],SiSinum)
						xSit=np.append(xSit[atomnum[0]*3:],atom[:atomnum[0]])
					
				
				linecount=1
				timecount+=1
			else:
				atom[linecount-1]=list(map(float,line.split()))
				linecount+=1
				continue
	return atomnum[0]
	
def plot(Sinum):
	#順番はSiO,SiSiとなっている
	start = time.time()
	with open("Sitimecord.dat","r") as Sc:
		print(Sc.readline())
		print(Sc.readline())
		time=0.001#ps
		count=0
		SiSi=np.array([])
		SiO=np.array([])
		while True:
			line=Sc.readline()
			if(line==""):
				print("OK")
				break
			temp=line.split()
			for i in temp:
				t1,t2=list(map(float,i.split(",")))
				SiSi=np.append(SiSi,int(t2))
				SiO=np.append(SiO,int(t1))
			count+=1
			time+=0.001
		SiSi=SiSi.reshape(int(len(SiSi)/Sinum),Sinum).T
		SiO=SiO.reshape(int(len(SiO)/Sinum),Sinum).T
		time=np.linspace(0.001,time,num=len(SiSi[0]))
		for i in range(len(SiSi)):
			plt.plot(time,SiSi[i])
			plt.savefig("SiSi"+str(i)+".png",dpi=150)
			plt.clf()
		for i in range(len(SiO)):
			plt.plot(time,SiO[i])
			plt.savefig("SiO"+str(i)+".png",dpi=150)
			plt.clf()
	elapsed_time = time.time() - start
	print ("elapsed_time:{0}[sec]".format(elapsed_time))
if __name__ =="__main__":
	if not os.path.exists("./change_cord"):
		os.mkdir("./change_cord")
	os.chdir("./change_cord")
	cord_time_change(3.0,2.0,1.6)
	plot(24)