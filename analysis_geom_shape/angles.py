from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math

def angles(traj):
    topology = traj.topology
    nres=traj.n_residues
    nch=traj.n_chains
    h,w = 501000, 250;
    ang0=np.zeros((h, w)) 
    ang1=np.zeros((h, w)) 
    ang2=np.zeros((h, w))
    ang3=np.zeros((h, w))
    ang4=np.zeros((h, w))
    ang5=np.zeros((h, w))
    ang6=np.zeros((h, w))
    ang7=np.zeros((h, w))
    
    k0=-1
    k1=-1
    k2=-1
    k3=-1
    k4=-1
    k5=-1
    k6=-1
    k7=-1
    
    for ii,rr in enumerate(topology.residues):
        
        bb1='%4d' % rr.index
        bb2='%4d' % rr.chain.index
        print(ii,rr)
        idxs = [[0,0,0] for x in range(10)]

            
        if(ii < nres-1):
            rr_p = topology.residue(ii+1)
            if(rr_p.index -1 == rr.index and rr_p.chain.index == rr.chain.index):
                try:
                    idxs[0]=([rr.atom("O3'").index,rr_p.atom("P").index, rr_p.atom("O2'").index])            
                    k0=k0+1
                    dd0=md.compute_angles(traj,[idxs[0]],periodic=False)
                    nsize=dd0.size
                    for ll in range(0,nsize):
                        ang0[ll][k0]=dd0[ll]
                except:
                    sys.exit("Error in the file to compute O3'  P O2' ang")

                    
                try:
                    idxs[1]=([rr.atom("O2'").index,rr_p.atom("P").index, rr_p.atom("O5'").index])            
                    k1=k1+1
                    dd1=md.compute_angles(traj,[idxs[1]],periodic=False)
                    nsize=dd1.size
                    for ll in range(0,nsize):
                        ang1[ll][k1]=dd1[ll]
                except:
                    sys.exit("Error in the file to compute O2'  P O5' ang")

                try:
                    idxs[2]=([rr.atom("O5'").index,rr_p.atom("P").index, rr_p.atom("O2'").index])            
                    k2=k2+1
                    dd2=md.compute_angles(traj,[idxs[2]],periodic=False)
                    nsize=dd2.size
                    for ll in range(0,nsize):
                        ang2[ll][k2]=dd2[ll]
                except:
                    sys.exit("Error in the file to compute  O5 P  O2 ang")

                try:
                    idxs[3]=([rr.atom("C1'").index,rr_p.atom("P").index, rr_p.atom("C1'").index])            
                    k3=k3+1
                    dd3=md.compute_angles(traj,[idxs[3]],periodic=False)
                    nsize=dd3.size
                    for ll in range(0,nsize):
                        ang3[ll][k3]=dd3[ll]
                except:
                    sys.exit("Error in the file to compute  C1  P  C2 ang")
                    
                try:
                    idxs[4]=([rr.atom("C2").index,rr_p.atom("P").index, rr_p.atom("C2").index])            
                    k4=k4+1
                    dd4=md.compute_angles(traj,[idxs[4]],periodic=False)
                    nsize=dd4.size
                    for ll in range(0,nsize):
                        ang4[ll][k4]=dd4[ll]
                except:
                    sys.exit("Error in the file to compute  C5 P  C2 ang")


                try:
                    idxs[5]=([rr.atom("C2").index,rr.atom("C1'").index, rr_p.atom("P").index])            
                    k5=k5+1
                    dd5=md.compute_angles(traj,[idxs[5]],periodic=False)
                    nsize=dd5.size
                    for ll in range(0,nsize):
                        ang5[ll][k5]=dd5[ll]
                except:
                    sys.exit("Error in the file to compute  C2 C1  P ang")
                    
                try:
                    idxs[6]=([rr.atom("C1'").index,rr.atom("C4'").index, rr_p.atom("P").index])            
                    k6=k6+1
                    dd6=md.compute_angles(traj,[idxs[6]],periodic=False)
                    nsize=dd6.size
                    for ll in range(0,nsize):
                        ang6[ll][k6]=dd6[ll]
                except:
                    sys.exit("Error in the file to compute  C1 C4  P ang")
                    
                try:
                    idxs[7]=([rr.atom("C2").index,rr.atom("C4'").index, rr_p.atom("P").index])            
                    k7=k7+1
                    dd7=md.compute_angles(traj,[idxs[7]],periodic=False)
                    nsize=dd7.size
                    for ll in range(0,nsize):
                        ang7[ll][k7]=dd7[ll]
                except:
                    sys.exit("Error in the file to compute  C2 C4  P ang")
                    
    conv=180.0/math.pi        
    np.savetxt('ang_O3_P_O2.dat',ang0[0:ll+1,0:k0+1]*conv,fmt='%4.2f')
    np.savetxt('ang_O2_P_O5.dat',ang1[0:ll+1,0:k1+1]*conv,fmt='%4.2f')
    np.savetxt('ang_O5_P_O2.dat',ang2[0:ll+1,0:k2+1]*conv,fmt='%4.2f')
    np.savetxt('ang_C1_P_C1.dat',ang3[0:ll+1,0:k3+1]*conv,fmt='%4.2f')
    np.savetxt('ang_C2_P_C2.dat',ang4[0:ll+1,0:k4+1]*conv,fmt='%4.2f')
    np.savetxt('ang_C2_C1_P.dat',ang5[0:ll+1,0:k5+1]*conv,fmt='%4.2f')
    np.savetxt('ang_C1_C4_P.dat',ang6[0:ll+1,0:k6+1]*conv,fmt='%4.2f')
    np.savetxt('ang_C2_C4_P.dat',ang7[0:ll+1,0:k7+1]*conv,fmt='%4.2f')
    
