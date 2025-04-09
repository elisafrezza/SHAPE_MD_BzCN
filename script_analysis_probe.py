import mdtraj as md
import numpy as np
import sys
import subprocess
import math

from findaxis import findaxis


typetop=input("Enter type of trajectory [pdb/xtc] ");
if(typetop=='pdb'):
    name=input("Enter name trajctory: ");
    traj1 = md.load(name+'.pdb')
elif(typetop=='xtc'):
    name=input("Enter name trajctory: ");
    topname=input("Enter name topology with extension: ");
    traj1 = md.load(name+'.xtc',top=topname)
print(traj1)
topology = traj1.topology
nres=traj1.n_residues
print(nres)
print(topology)
nsnap=traj1.n_frames
print(nsnap)
#traj=traj1
traj=traj1.image_molecules()
n_nucl=int(input("Provide the number of nucleotide you want to investigate for the nucleobase "));
n_nucl1=int(input("Provide the number of nucleotide studied in the Umbrella Sampling "));
probe=str(input("Provide the name of the probe (1M7/BzCN/NMIA)"));
###COM and vectors PROBE
xcom=np.zeros(nsnap)
ycom=np.zeros(nsnap)
zcom=np.zeros(nsnap)
vec1=np.zeros([nsnap,3])
vec2=np.zeros([nsnap,3])
normp=np.zeros([nsnap,3])
vec1l=np.zeros(3)
vec2l=np.zeros(3)

###COM and vectors PROBE
xcombase=np.zeros(nsnap)
ycombase=np.zeros(nsnap)
zcombase=np.zeros(nsnap)
vec1base=np.zeros([nsnap,3])
vec2base=np.zeros([nsnap,3])
normpbase=np.zeros([nsnap,3])
unit2pr=np.zeros([nsnap,3])
##PROBE
if(probe=='1M7'):
    atomid=topology.select("(resid 14 and name C1) or (resid 14 and name O1) or (resid 14 and name C2)  or (resid 14 and name C3) or (resid 14 and name C5) or (resid 14 and name C6)  or (resid 14 and name C7) or (resid 14 and name C8) or (resid 14 and name C4) or (resid 14 and name N1)"  )
    #print(atomid)
    #print(len(atomid))
    atomP1=topology.select(" (resid 14 and name N1)"  )
    atomP2=topology.select(" (resid 14 and name O1)"  )
    # C1,O1,C2,C3,C5,C6,C7,C8, C4, N1
    #print(traj.xyz[:,atomid,:])
    natp=len(atomid)
elif(probe=='BzCN'):
    atomid=topology.select("(resid 14 and name C1) or (resid 14 and name C2)  or (resid 14 and name C3) or (resid 14 and name C4) or (resid 14 and name C5)  or (resid 14 and name C6) "  )
    #print(atomid)
    #print(len(atomid))
    atomP1=topology.select(" (resid 14 and name C5)"  )
    atomP2=topology.select(" (resid 14 and name C6)"  )
    natp=len(atomid)

elif(probe=='NMIA'):
    atomid=topology.select("(resid 14 and name C1) or (resid 14 and name O2) or (resid 14 and name C5)  or (resid 14 and name N1) or (resid 14 and name C3) or (resid 14 and name C6)  or (resid 14 and name C8) or (resid 14 and name C7) or (resid 14 and name C4) or (resid 14 and name C2)"  )
    #print(atomid)
    #print(len(atomid))
    
    atomP1=topology.select(" (resid 14 and name N1)"  )
    atomP2=topology.select(" (resid 14 and name C1)"  )
    # C1,O1,C2,C3,C5,C6,C7,C8, C4, N1
    #print(traj.xyz[:,atomid,:])
    natp=len(atomid)
    #print('natp=',natp)
    # C1,O1,C2,C3,C5,C6,C7,C8, C4, N1
    
## Determination COM and  axis Probe
normloc=np.zeros(3)
for i in range(0,nsnap):
    
    xcom[i]=np.sum(traj.xyz[i,atomid,0])/natp
    ycom[i]=np.sum(traj.xyz[i,atomid,1])/natp
    zcom[i]=np.sum(traj.xyz[i,atomid,2])/natp

    vec1l[0]=traj.xyz[i,atomP1,0]-xcom[i]
    vec2l[0]=traj.xyz[i,atomP2,0]-xcom[i]
    vec1l[1]=traj.xyz[i,atomP1,1]-ycom[i]
    vec2l[1]=traj.xyz[i,atomP2,1]-ycom[i]
    vec1l[2]=traj.xyz[i,atomP1,2]-zcom[i]
    vec2l[2]=traj.xyz[i,atomP2,2]-zcom[i]
    #print(np.sqrt(vec1l[0]**2+vec1l[1]**2+vec1l[2]**2))
    vec1[i,0]=vec1l[0]/np.linalg.norm(vec1l)
    vec1[i,1]=vec1l[1]/np.linalg.norm(vec1l)
    vec1[i,2]=vec1l[2]/np.linalg.norm(vec1l)
    #print(np.sqrt(vec1[i,0]**2+vec1[i,1]**2+vec1[i,2]**2))
    vec2[i,0]=vec2l[0]/np.linalg.norm(vec2l)
    vec2[i,1]=vec2l[1]/np.linalg.norm(vec2l)
    vec2[i,2]=vec2l[2]/np.linalg.norm(vec2l)
## norm vector
    normloc=np.cross(vec1[i,:],vec2[i,:])
    normp[i,0]=normloc[0]/np.linalg.norm(normloc)
    normp[i,1]=normloc[1]/np.linalg.norm(normloc)
    normp[i,2]=normloc[2]/np.linalg.norm(normloc)
    unit2pr[i,:]=np.cross(vec1[i,:],normp[i,:])

########################
####Nucleotide
## Nucleobase
# Nucelotide under investigation
nres_loc=n_nucl-56
res=str(topology.residue(nres_loc))[0:1]
if(res=='G' or res=='A'):

    atomidbase=topology.select("(resid "+str(nres_loc)+" and name N9) or (resid "+str(nres_loc)+" and name C4 ) or (resid "+str(nres_loc)+" and name N3 ) or (resid "+str(nres_loc)+" and name C2 ) or (resid "+str(nres_loc)+" and name N1 ) or (resid "+str(nres_loc)+" and name C6 ) or (resid "+str(nres_loc)+" and name C5 ) or (resid "+str(nres_loc)+" and name N7 ) or (resid "+str(nres_loc)+" and name C8)")

    natbase=len(atomidbase)
    stringa='resid ' 
    atom1=topology.select("(resid "+str(nres_loc)+" and name N3)")
    #print(atom1)
    atom2=topology.select("(resid "+str(nres_loc)+" and name C6)")


unit2=np.zeros([nsnap,3])
norml=np.zeros(3)
displbase=np.zeros([nsnap,3])
#-----------
## Determination COM and  axis 
for i in range(0,nsnap):
    
    xcombase[i]=np.sum(traj.xyz[i,atomidbase,0])/natbase
    ycombase[i]=np.sum(traj.xyz[i,atomidbase,1])/natbase
    zcombase[i]=np.sum(traj.xyz[i,atomidbase,2])/natbase

    #print(traj.xyz[i,atom1,0]-xcombase[i])

    vec1l[0]=traj.xyz[i,atom1,0]-xcombase[i]
    vec1l[1]=traj.xyz[i,atom1,1]-ycombase[i]
    vec1l[2]=traj.xyz[i,atom1,2]-zcombase[i]
    
    vec2l[0]=traj.xyz[i,atom2,0]-xcombase[i]
    vec2l[1]=traj.xyz[i,atom2,1]-ycombase[i]
    vec2l[2]=traj.xyz[i,atom2,2]-zcombase[i]

    #print(np.sqrt(vec1l[0]**2+vec1l[1]**2+vec1l[2]**2))

    vec1base[i,0]=vec1l[0]/np.linalg.norm(vec1l)
    vec1base[i,1]=vec1l[1]/np.linalg.norm(vec1l)
    vec1base[i,2]=vec1l[2]/np.linalg.norm(vec1l)
    
    #print(np.sqrt(vec1[i,0]**2+vec1[i,1]**2+vec1[i,2]**2))
    vec2base[i,0]=vec2l[0]/np.linalg.norm(vec2l)
    vec2base[i,1]=vec2l[1]/np.linalg.norm(vec2l)
    vec2base[i,2]=vec2l[2]/np.linalg.norm(vec2l)
    
    norml=np.cross(vec1base[i,:],vec2base[i,:])
    normpbase[i,0]=norml[0]/np.linalg.norm(norml)
    normpbase[i,1]=norml[1]/np.linalg.norm(norml)
    normpbase[i,2]=norml[2]/np.linalg.norm(norml)

    unit2[i,:]=np.cross(vec1base[i,:],normpbase[i,:])




###Allocation vectors
angaxispb=np.zeros([nsnap,4])
distcompb=np.zeros(nsnap)
stackbase=np.zeros(nsnap)
dist=np.zeros(3)
dr1=np.zeros(3)
ww=np.zeros(3)
ua=np.zeros(3)
matrixq=np.zeros([3,3])
b1=np.zeros([3,3])
b2=np.zeros([3,3])
theta=np.zeros(nsnap)
ang=np.zeros([nsnap,3])
for i in range(0,nsnap):
    angaxispb[i,0]=math.degrees(np.arccos(np.absolute(np.dot(normp[i,:],normpbase[i,:]))))
    dx=xcom[i]-xcombase[i]
    dy=ycom[i]-ycombase[i]
    dz=zcom[i]-zcombase[i]
    dr1[0]=dx
    dr1[1]=dy
    dr1[2]=dz
    distcompb[i]=math.sqrt(dx**2+dy**2+dz**2)

    dist[0]=dx/distcompb[i]
    dist[1]=dy/distcompb[i]
    dist[2]=dz/distcompb[i]
    #print(math.sqrt(dist[0]**2+dist[1]**2+dist[2]**2))
    #print(np.dot(vec1[i,:],vec1base[i,:]))
    #angle between the plates
    angaxispb[i,1]=math.degrees(np.arccos(np.dot(vec1[i,:],vec1base[i,:])))
    #relative angle for the displacement 
    angaxispb[i,2]=math.degrees(np.arccos(np.dot(vec1[i,:],dist)))
    angaxispb[i,3]=math.degrees(np.arccos(np.dot(vec1base[i,:],dist)))
    #print(angaxispb[i,:])
    if(distcompb[i]<0.5):
        if(angaxispb[i,0]<45.0):
            stackbase[i]=1.0
    displbase[i,0]=np.dot(vec1base[i,:],dr1)*10
    displbase[i,1]=np.dot(unit2[i,:],dr1)*10
    displbase[i,2]=np.dot(normpbase[i,:],dr1)*10

    #Q matrix
    b1[0,:]=vec1base[i,:]
    b1[1,:]=unit2[i,:]
    b1[2,:]=normpbase[i,:]
    b2[0,:]=vec1[i,:]
    b2[1,:]=unit2pr[i,:]
    b2[2,:]=normp[i,:]

    matrixq=np.matmul(np.matrix.transpose(b1),b2)
    #I_p,P=np.linalg.eig(matrixq)
    #print(I_p)
    #print(P)
    ct=(np.trace(matrixq)-1)/2
    if np.abs(ct) > 1:
        ct = np.sign(ct)
    theta=np.arccos(ct)
    if(theta< 1e-6):
        print('WARNING')
    elif(math.pi-theta< 1e-6):
        aaa=matrixq+1.0
        tol1 = 1e-10
        tol2 = 1e-6
        max_iter = 40

        # Initialize gamma, usc, and it
        gamma = 0.0
        ww  = 0.0
        it = -1
        gamma, ww, it = findaxis(aaa, tol1, tol2, max_iter, gamma, ww, it)
    else:
        ww[0]=matrixq[1,2]-matrixq[2,1]
        ww[1]=matrixq[2,0]-matrixq[0,2]
        ww[2]=matrixq[0,1]-matrixq[1,0]
    ua[0]=ww[0]/np.linalg.norm(ww)
    ua[1]=ww[1]/np.linalg.norm(ww)
    ua[2]=ww[2]/np.linalg.norm(ww)
    thetadeg=theta*180/math.pi
    ag1=np.dot(ua,vec1base[i,:])*thetadeg
    ag2=np.dot(ua,unit2[i,:])*thetadeg
    ag3=np.dot(ua,normpbase[i,:])*thetadeg
    ang[i,0]=np.dot(ua,b1[0,:])*thetadeg
    ang[i,1]=np.dot(ua,b1[1,:])*thetadeg
    ang[i,2]=np.dot(ua,b1[2,:])*thetadeg


np.savetxt('angle_probe_'+probe+'_base_'+str(n_nucl)+'.dat',angaxispb,fmt='%6.2f')
np.savetxt('distance_com_probe_'+probe+'_base_'+str(n_nucl)+'.dat',distcompb*10.0,fmt='%5.2f')
np.savetxt('stack_'+probe+'_base_'+str(n_nucl)+'.dat',stackbase,fmt='%5.2f')
np.savetxt('displ_com_probe_'+probe+'_base_'+str(n_nucl)+'.dat',displbase,fmt='%5.2f')

################
##Ribose
#print(atomidbase)
nres_loc=n_nucl-56

at=topology.select
#index

atomidrib=[]
for i in range(0,5):
    atomidrib.append(0)
    
atomidrib[0]=topology.residue(nres_loc).atom("O4'").index
atomidrib[1]=topology.residue(nres_loc).atom("C4'").index
atomidrib[2]=topology.residue(nres_loc).atom("C3'").index
atomidrib[3]=topology.residue(nres_loc).atom("C2'").index
atomidrib[4]=topology.residue(nres_loc).atom("C1'").index

#print(atomidrib)
atomidrib1=topology.select("(resid "+str(nres_loc)+" and index "+str(atomidrib[0])+") or (resid "+str(nres_loc)+" and index "+str(atomidrib[1])+") or (resid "+str(nres_loc)+" and index "+str(atomidrib[2])+") or (resid "+str(nres_loc)+" and index "+str(atomidrib[3])+") or (resid "+str(nres_loc)+" and index "+str(atomidrib[4])+") ")
#atomidrib=topology.select("(resid "+str(nres_loc)+" and name C1')")

atom1r=topology.select("(resid "+str(nres_loc)+" and index "+str(atomidrib[0])+")")

atom2r=topology.select("(resid "+str(nres_loc)+" and index "+str(atomidrib[1])+")")


massesrib=[14,14,14,14,16]

###Allocation

#angaxispb=np.zeros([nsnap,3])
#distcompb=np.zeros(nsnap)
comrib=np.zeros(3)
unitvec=np.zeros(3)
angaxispr=np.zeros([nsnap,4])
distcompr=np.zeros(nsnap)
vec1rib=np.zeros(3)
stackrib=np.zeros(nsnap)
vec2rib=np.zeros(3)
dr2=np.zeros(3)
displrib=np.zeros([nsnap,3])
angr=np.zeros([nsnap,3])
aaa=np.zeros([3,3])
test=np.zeros(3)
vec2l=np.zeros(3)
vec2t=np.zeros(3)
testl=np.zeros(3)
for i in range(0,nsnap):
    comrib[0]=np.sum(traj.xyz[i,atomidrib1,0])/5.0
    comrib[1]=np.sum(traj.xyz[i,atomidrib1,1])/5.0
    comrib[2]=np.sum(traj.xyz[i,atomidrib1,2])/5.0
   # print('COM RIBOSE',comrib*10,traj.xyz[i,atomidrib1[0],0])
    newcoor=traj.xyz[i,atomidrib1,:]-comrib
    ref2=traj.xyz[i,1,:]-comrib
    #print(newcoor)
    #traj.xyz=newcoor
    #printtraj.xyz)
    x=newcoor[:,0]
    y=newcoor[:,1]
    z=newcoor[:,2]
    masses=1.0
    Ixx = np.sum(masses * (y**2 + z**2))
    Iyy = np.sum(masses * (x**2 + z**2))
    Izz = np.sum(masses * (x**2 + y**2))
    Ixy = -np.sum(masses * x * y)
    Iyz = -np.sum(masses * y * z)
    Ixz = -np.sum(masses * x * z)
    II = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    I_p,P=np.linalg.eig(II)

    #print(I_p,P)
    #print(I_p[0])
    aa=1000000
    ind=0
    for i2  in range(0,3):
        if(I_p[i2]<aa):
            ind=i2
            aa=I_p[i2]
            #print(aa,ind)
    axes=P[:,0]

   
    #print('test',np.dot(diff/np.linalg.norm(diff),ref2/np.linalg.norm(ref2)),'axes',axes*10+comrib*10)
  
    #print(axes,comrib,(unitvec+comrib)*10,comrib*10)
    angaxispr[i,0]=math.degrees(np.arccos(np.absolute(np.dot(normp[i,:],unitvec))))

    dx=xcom[i]-comrib[0]
    dy=ycom[i]-comrib[1]
    dz=zcom[i]-comrib[2]
    dr2[0]=dx
    dr2[1]=dy
    dr2[2]=dz
    #print(dr2)
    distcompr[i]=math.sqrt(dx**2+dy**2+dz**2)
    dist[0]=dx/np.linalg.norm(distcompr[i])#distcompr[i]
    dist[1]=dy/np.linalg.norm(distcompr[i])#distcompr[i]
    dist[2]=dz/np.linalg.norm(distcompr[i])#distcompr[i]
    #print(distcompr[i],np.linalg.norm(distcompr[i]),np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2))
    #print(dist)
    #angle between the plates
    vec1l[0]=traj.xyz[i,atom1r,0]-comrib[0]
    vec1l[1]=traj.xyz[i,atom1r,1]-comrib[1]
    vec1l[2]=traj.xyz[i,atom1r,2]-comrib[2]
    vec1rib[0]=vec1l[0]/np.linalg.norm(vec1l)
    vec1rib[1]=vec1l[1]/np.linalg.norm(vec1l)
    vec1rib[2]=vec1l[2]/np.linalg.norm(vec1l)
    vec2l[0]=traj.xyz[i,atom2r,0]-comrib[0]
    vec2l[1]=traj.xyz[i,atom2r,1]-comrib[1]
    vec2l[2]=traj.xyz[i,atom2r,2]-comrib[2]
    vec2t[0]=vec2l[0]/np.linalg.norm(vec2l)
    vec2t[1]=vec2l[1]/np.linalg.norm(vec2l)
    vec2t[2]=vec2l[2]/np.linalg.norm(vec2l)
    testl=np.cross(vec1rib,vec2t)
    #/np.linalg.norm((np.cross(vec1rib,vec2rib)))
    test[0]=testl[0]/np.linalg.norm(testl)
    test[1]=testl[1]/np.linalg.norm(testl)
    test[2]=testl[2]/np.linalg.norm(testl)
    diff=axes-comrib
   # print(test,axes,np.dot(axes,test))
    if(np.dot(axes,test)>0):
#    if(np.dot(diff/np.linalg.norm(diff),ref2/np.linalg.norm(ref2))<0):
        axes=-axes
    unitvec=axes/np.linalg.norm(axes)
    #3D perpendicular axis
    vec2rib=np.cross(vec1rib,unitvec)
    

    angaxispr[i,1]=math.degrees(np.arccos(np.dot(vec1[i,:],vec1rib)))

    angaxispr[i,2]=math.degrees(np.arccos(np.dot(vec1[i,:],dist)))
    #print(np.dot(vec1rib,dist))
    angaxispr[i,3]=math.degrees(np.arccos(np.dot(vec1rib,dist)))

    displrib[i,0]=np.dot(vec1rib,dr2)*10
    displrib[i,1]=np.dot(vec2rib,dr2)*10
    displrib[i,2]=np.dot(unitvec,dr2)*10
    #print(angaxispr)
    if(distcompr[i]<0.5):
        if(angaxispr[i,0]<45.0):
            stackrib[i]=1.0
            
        #Q matrix
    b1[0,:]=vec1rib
    b1[1,:]=vec2rib
    b1[2,:]=unitvec
    b2[0,:]=vec1[i,:]
    b2[1,:]=unit2pr[i,:]
    b2[2,:]=normp[i,:]

    matrixq=np.matmul(np.matrix.transpose(b1),b2)
    
    #I_p,P=np.linalg.eig(matrixq)
    #print(I_p)
    #print(P)
    ct=(np.trace(matrixq)-1)/2
    #print(ct)
    if np.abs(ct) > 1:
        ct = np.sign(ct)
    theta=np.arccos(ct)
    #print(theta,ct)
    thetadeg=theta*180/math.pi
    
    if(theta< 1e-6):
          print('WARNING')
    elif(math.pi-theta< 1e-6):
        aaa=matrixq+1.0
        tol1 = 1e-10
        tol2 = 1e-6
        max_iter = 40

        # Initialize gamma, usc, and it
        gamma = 0.0
        ww  = 0.0
        it = -1
        gamma, ww, it = findaxis(aaa, tol1, tol2, max_iter, gamma, ww, it)
    else:
        ww[0]=matrixq[1,2]-matrixq[2,1]
        ww[1]=matrixq[2,0]-matrixq[0,2]
        ww[2]=matrixq[0,1]-matrixq[1,0]
        
    ua[0]=ww[0]/np.linalg.norm(ww)
    ua[1]=ww[1]/np.linalg.norm(ww)
    ua[2]=ww[2]/np.linalg.norm(ww) 
    ag1=np.dot(ua,vec1base[i,:])*thetadeg
    ag2=np.dot(ua,unit2[i,:])*thetadeg
    ag3=np.dot(ua,normpbase[i,:])*thetadeg
    angr[i,0]=np.dot(ua,b1[0,:])*thetadeg
    angr[i,1]=np.dot(ua,b1[1,:])*thetadeg
    angr[i,2]=np.dot(ua,b1[2,:])*thetadeg
    
#print(displbase,displrib,dr2,vec1rib)
np.savetxt('angle_probe_'+probe+'_ribose_'+str(n_nucl)+'.dat',angaxispr,fmt='%6.2f')
np.savetxt('distance_com_probe_'+probe+'_ribose_'+str(n_nucl)+'.dat',distcompr*10.0,fmt='%5.2f')
np.savetxt('displ_com_probe_'+probe+'_ribose_'+str(n_nucl)+'.dat',displrib,fmt='%5.2f')
np.savetxt('stack_'+probe+'_ribose_'+str(n_nucl)+'.dat',stackrib,fmt='%5.2f')

