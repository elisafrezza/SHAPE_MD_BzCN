from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math

def puckering(traj):
    topology = traj.topology
    nres=traj.n_residues
    h,w = 1000001, 5;
    dih0=np.zeros((h, w)) # O5  C5  CA  CY
    puck=np.zeros((h, 200))
    pha2=np.zeros((h, 200))
    crd=57.29577951
    k0=-1
    k1=-1
    k2=-1
    k3=-1
    k4=-1
    nchref=-1
    kref=-1
    sugt=[2,8,4,5,1,7,3,9,0,6]
    for ii,rr in enumerate(topology.residues):
        
        bb1='%4d' % rr.index
        bb2='%4d' % rr.chain.index
        #print(ii,rr)
        idxs = [[0,0,0,0] for x in range(15)]
        if(nchref !=rr.chain.index):
            nchref=rr.chain.index
            kref=0
        else:
            kref=kref+1

        try:
            idxs[0] =([rr.atom("C1'").index, rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index])
            k0=k0+1
        except:
            sys.exit("Error in the file dihedral v1")  
      
        
        ang=md.compute_dihedrals(traj,[idxs[0]],periodic=False)
        #print(ang)
        #print(ang[1])
        nsize=ang.size
    
        for ll in range(0,nsize):
            dih0[ll][0]=ang[ll]
       
        try:
            idxs[1] =([rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index])
            k1=k1+1
        except:
            sys.exit("Error in the file dihedral v2 ")
            
                
        ang1=md.compute_dihedrals(traj,[idxs[1]],periodic=False)
        #        print(ang)
        #print(ang[1])
        nsize1=ang1.size
        
        for ll in range(0,nsize1):
            dih0[ll][1]=ang1[ll]

        try:
            idxs[2] =([rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index])
            k2=k2+1
        except:
            sys.exit("Error in the file dihedral v3")
            
    
                
        ang2=md.compute_dihedrals(traj,[idxs[2]],periodic=False)
        #        print(ang)
        #print(ang[1])
        nsize2=ang2.size
    
        for ll in range(0,nsize2):
            dih0[ll][2]=ang2[ll]


            
        try:
            idxs[3] =([rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index,rr.atom("C2'").index])
            k3=k3+1
        except:
            sys.exit("Error in the file dihedral v4")

        ang3=md.compute_dihedrals(traj,[idxs[3]],periodic=False)
        #        print(ang)
        #print(ang3[0])
        nsize3=ang3.size
    
        for ll in range(0,nsize3):
            dih0[ll][3]=ang3[ll]

        try:
            idxs[4] =([rr.atom("O4'").index, rr.atom("C1'").index,rr.atom("C2'").index,rr.atom("C3'").index])
            k4=k4+1
        except:
            sys.exit("Error in the file dihedral v4")

        ang4=md.compute_dihedrals(traj,[idxs[4]],periodic=False)
        #        print(ang)
        #print(ang3[0])
        nsize4=ang4.size
    
        for ll in range(0,nsize4):
            
            dih0[ll][4]=ang4[ll]

        
        for ll  in range(0,nsize4):
            a=0
            b=0
            for tt in range(0,5):
               a=a+dih0[ll,tt]*math.cos(144*tt/crd)*crd
               b=b-dih0[ll,tt]*math.sin(144*tt/crd)*crd
        
            a=a*0.4
            b=b*0.4
            
            

            amp=math.sqrt(a*a+b*b)
            
            
            if(amp>0):
                cp=a/amp
                sp=b/amp
             
                if(abs(cp)>1):
                    cp=copysign(1,cp)
                pha=math.acos(cp)*crd
               
                if(sp<0):
                    pha=-pha
            else:
                pha=0.0
         
            if(pha < 0):
                pha=pha+360.0
         
            puck[ll,k0]=sugt[int(pha/36)]
            print(puck[ll,k0],pha)
            pha2[ll,k0]=pha
             
        conv=180.0/math.pi

        
    np.savetxt('puck_trj.dat',puck[0:ll+1,0:k0+1],fmt='%4d')
    np.savetxt('phase_trj.dat',pha2[0:ll+1,0:k0+1],fmt='%4d')
        #np.savetxt('dih_NbC1pC3pO2p_trj.dat',dih1[0:ll+1,0:k1+1]*conv,fmt='%4.2f')
        #np.savetxt('dih_NbC1pC4pO2p_trj.dat',dih2[0:ll+1,0:k2+1]*conv,fmt='%4.2f')
