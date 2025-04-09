import mdtraj as md
import numpy as np
import sys
import subprocess
import math
from datetime import datetime


typetop=input("Enter type of trajectory [pdb/xtc] ");
if(typetop=='pdb'):
    name=input("Enter name trajctory: ");
    traj = md.load(name+'.pdb')
elif(typetop=='xtc'):
    name=input("Enter name trajctory: ");
    topname=input("Enter name topology with extension: ");
    traj = md.load(name+'.xtc',top=topname)


start = datetime.now()

print(traj)
topology = traj.topology
nres=traj.n_residues
print(nres)
print(topology)
nsnap=traj.n_frames
print(nsnap)



CUTOFF=1.0

index=topology.select("resid 14 ")
nprobe=len(index) #topology.n_atoms #475
istart=index[0]
iend=index[nprobe-1]+1

##Allocation
CV=np.zeros([nsnap,nprobe])
CVw=np.zeros([nsnap,nprobe])

iloc=-1
for i in range(istart,iend):
    iloc=iloc+1
    print('atom ',i)
    atomO2=i
    #id_O2=topology.select("(resid "+str(i)+" and index "+str(atomO2)+")")
    id_O2=topology.select("(index "+str(atomO2)+")")
    #    id_all=topology.select("not protein and not (resid "+str(i)+" and index "+str(atomO2)+")")
    id_all=topology.select(" not (index "+str(atomO2)+") or (not resid 14)")
    #print(id_O2)
    #print(id_all)
    pairs = np.array([(i,j) for  i in id_O2  for j in id_all])
    #print(pairs)
    pairs_distances_traj = md.compute_distances(traj,pairs)
    #print(pairs_distances_traj)
    for ii in range(0,traj.n_frames):
        
        contacts = pairs[pairs_distances_traj[ii,:] < CUTOFF]
        index_cut=np.unique(contacts[:,1])
        
        den=0.0
        dxx=0.0
        dyy=0.0
        dzz=0.0
        dxx2=0.0
        dyy2=0.0
        dzz2=0.0
        ###CV unit vector
        axx=traj.xyz[ii,i,0]-traj.xyz[ii,index_cut,0]
        ayy=traj.xyz[ii,i,1]-traj.xyz[ii,index_cut,1]
        azz=traj.xyz[ii,i,2]-traj.xyz[ii,index_cut,2]
        dddx=axx/np.sqrt(axx**2+ayy**2+azz**2)
        dddy=ayy/np.sqrt(axx**2+ayy**2+azz**2)
        dddz=azz/np.sqrt(axx**2+ayy**2+azz**2)
        dd4=np.sqrt(np.sum(dddx)**2+np.sum(dddy)**2+np.sum(dddz)**2)
        CV[ii,iloc]=1-dd4/len(index_cut)
        ### CVw
        num1=np.sqrt(np.sum(np.sum(axx))**2+np.sum(np.sum(ayy))**2+np.sum(np.sum(azz)**2))
        den1=np.sum(np.sqrt(axx**2+ayy**2+azz**2))
        CVw[ii,iloc]=1-num1/den1



        
np.savetxt('CV.dat',CV,fmt='%6.3f')
np.savetxt('CVw.dat',CVw,fmt='%6.3f')

CVstat=np.zeros([nprobe,2])
CVwstat=np.zeros([nprobe,2])

for i in range(0,nprobe):
    CVstat[i,0]=np.average(CV[:,i])
    CVstat[i,1]=np.std(CV[:,i])
    CVwstat[i,0]=np.average(CVw[:,i])
    CVwstat[i,1]=np.std(CVw[:,i])


np.savetxt('CV_stat.dat',CVstat,fmt='%6.3f')
np.savetxt('CVw_stat.dat',CVwstat,fmt='%6.3f')

difference = datetime.now() - start
print(difference)
