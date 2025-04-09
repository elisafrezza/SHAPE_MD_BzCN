from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math

#name=input("Enter name trajctory")
#file1=open(name+'.pdb','r')
#file2=open('snap1.pdb','w')
#file4=open('cc.pdb','w')
#k=0
#k1=-1
#for line in file1:
#    if(line[0:4]=='ATOM'):
#        #print line
#        file2.write(line)
#    elif(line[0:6]=='ENDMDL'):
#        file2.close()
#        subprocess.call('/mnt/d11266a1-4030-4f79-933c-eabe2431b041/NucleicAcids/HiRE/RunStartUp/HiRE_fa2cg/fa2cg/FA2CG < snap1.pdb > snap1_CG.pdb',shell='True')
#        file3=open('snap1_CG.pdb','r')
#        for line2 in file3:
#            file4.write(line2)
#        file3.close()
#        file2=open('snap1.pdb','w')
user_input = input("Enter the analysis you want to performe in one line using commas -- To help, type help ")
 
if(user_input == 'help'):
    print('Possible analysis:')
    print('TORSIONS: Torsional analysis')
    print('BONDS: Bond analysis')
    print('ANGLES: Angle analysis')
    print('PUCKERING: Angle analysis')
    
    sys.exit()


input_list = user_input.split(',')
print(input_list)
nsizel=len(input_list)
print(nsizel)

typetop=input("Enter type of trajectory [pdb/xtc] ");
if(typetop=='pdb'):
    name=input("Enter name trajctory: ");
    traj = md.load(name+'.pdb')
elif(typetop=='xtc'):
    name=input("Enter name trajctory: ");
    topname=input("Enter name topology with extension: ");
    traj = md.load(name+'.xtc',top=topname)
print(traj)
topology = traj.topology
nres=traj.n_residues
print(nres)
print(topology)



for il in range(0,nsizel):
    if (input_list[il]=='TORSIONS'):
        import torsions
        torsions.torsions(traj)
        
    elif(input_list[il]=='BONDS'):
        
        import bonds
        bonds.bonds(traj)

    elif(input_list[il]=='ANGLES'):
        
        import angles
        angles.angles(traj)
    elif(input_list[il]=='PUCKERING'):
        
        import puckering
        puckering.puckering(traj)
