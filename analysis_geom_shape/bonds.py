from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math

def bonds(traj):
    topology = traj.topology
    nres=traj.n_residues
    nch=traj.n_chains
    h,w = 501000, 250;
    bond0=np.zeros((h, w)) 
    bond1=np.zeros((h, w)) 
    bond2=np.zeros((h, w))
    bond3=np.zeros((h, w))
    bond4=np.zeros((h, w))
    bond5=np.zeros((h, w))

    ntframes=traj.n_frames
    k0=-1
    k1=-1
    k2=-1
    k3=-1
    k4=-1
    k5=-1
    
    for ii,rr in enumerate(topology.residues):
        
        bb1='%4d' % rr.index
        bb2='%4d' % rr.chain.index
       # print(ii,rr)
        idxs = [[0,0] for x in range(6)]
 
        if(ii < nres-1):
            rr_p = topology.residue(ii+1)
            if(rr_p.index -1 == rr.index and rr_p.chain.index == rr.chain.index):
                try:
                    idxs[0]=([rr.atom("C2").index,rr_p.atom("C2").index])            
                    k0=k0+1
                    #print(k0,idxs[0])
                    dd0=md.compute_distances(traj,[idxs[0]],periodic=False)
                    nsize=dd0.size
                    #print(dd0)
                    for ll in range(0,nsize):
                        bond0[ll][k0]=dd0[ll]
                except:
                    sys.exit("Error in the file to compute C2_i C2_i+1  bond")
                try:
                    idxs[1]=([rr.atom("C1'").index, rr_p.atom("C1'").index])            
                    k1=k1+1
                    dd1=md.compute_distances(traj,[idxs[1]],periodic=False)
                    nsize=dd1.size
                    for ll in range(0,nsize):
                        bond1[ll][k1]=dd1[ll]
                except:
                    sys.exit("Error in the file to compute C1' C2'  bond")
                if(ii > 0 ):
                    try:
                        idxs[2]=([rr.atom("P").index, rr_p.atom("P").index])            
                        k2=k2+1
                        dd2=md.compute_distances(traj,[idxs[2]],periodic=False)
                        nsize=dd2.size
                        for ll in range(0,nsize):
                            bond2[ll][k2]=dd2[ll]
                    except:
                        sys.exit("Error in the file to compute P P  bond")
        if( ii < nres-1):
            rr_p = topology.residue(ii+1)
            indo1=rr_p.atom("OP1").index
            indo2=rr_p.atom("OP2").index
            indP=rr_p.atom("P").index
            indO2=rr.atom("O2'").index
            #print(indo1,indo2,indP,indO2)
           # print(traj.xyz[ifr,indo1,:])
            
            k3=k3+1
            k5=k5+1
            for ifr in range(0,ntframes):
                #dist1=traj.xyz[ifr,indo1,:]
                dr1=(traj.xyz[ifr,indo1,:]+traj.xyz[ifr,indo2,:]-traj.xyz[ifr,indP,:]-traj.xyz[ifr,indO2,:])
                dist1=np.sqrt(np.sum(dr1**2))
                #print(traj.xyz[:,indo1,0],traj.xyz[:,indo1,1],traj.xyz[:,indo1,2])
                #print(traj.xyz[ifr,indo1,:]+traj.xyz[ifr,indo2,:]-traj.xyz[ifr,indP,:])
                #print('dist',ifr,ii,dr1*10,dist1)
                bond3[ifr][k3]=dist1
    
                dr2a=(traj.xyz[ifr,indo1,:]-traj.xyz[ifr,indO2,:])
                dist2a=np.sqrt(np.sum(dr2a**2))
                dr2b=(traj.xyz[ifr,indo2,:]-traj.xyz[ifr,indO2,:])
                dist2b=np.sqrt(np.sum(dr2b**2))
                if(dist2a < dist2b):
                    bond5[ifr][k5]=dist2a
                else:
                    bond5[ifr][k5]=dist2b
        #print(rr..name)
        if(rr.name=='G' or rr.name=='A'):
            try:
                idxs[3]=([rr.atom("O2'").index,rr.atom("N3").index])            
                k4=k4+1
                #print(k4,idxs[0])
                dd4=md.compute_distances(traj,[idxs[3]],periodic=False)
                nsize=dd4.size
                #print(dd4)
                for ll in range(0,nsize):
                    bond4[ll][k4]=dd4[ll]
            except:
                sys.exit("Error in the file to compute O2'_i N3_i  bond")       
        else:
            try:
                idxs[3]=([rr.atom("O2'").index,rr.atom("O2").index])            
                k4=k4+1
                #print(k4,idxs[0])
                dd4=md.compute_distances(traj,[idxs[3]],periodic=False)
                nsize=dd4.size
                #print(dd4)
                for ll in range(0,nsize):
                    bond4[ll][k4]=dd4[ll]
            except:
                sys.exit("Error in the file to compute O2'_i O2_i  bond")   
        
    np.savetxt('dist_C2.dat',bond0[0:ll+1,0:k0+1]*10.0,fmt='%4.2f')
    np.savetxt('dist_C1.dat',bond1[0:ll+1,0:k1+1]*10,fmt='%4.2f')
    np.savetxt('dist_P.dat',bond2[0:ll+1,0:k2+1]*10,fmt='%4.2f')
    np.savetxt('dist_OP_O2p_bis.dat',bond3[0:ll+1,0:k3+1]*10,fmt='%4.2f')
    np.savetxt('dist_O2p_base.dat',bond4[0:ll+1,0:k4+1]*10,fmt='%4.2f')
    np.savetxt('dist_OP_O2p_min.dat',bond5[0:ll+1,0:k5+1]*10,fmt='%4.2f')
