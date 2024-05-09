#!/usr/bin/env python
#(c) Yuan Zhang in Department of Statistic, FSU
# Nov.1, 2017
# This script is used to get the key atoms on protein surface

import sys
import numpy as np
import pylab as py
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift


def read_file(file_name):       # skip alternate position version
    data = []
    count = -1
    residue = 0     # label to track if the residue changed
    resname = []
    with open (file_name) as f:
        for line in f:
            a = line
            # x [[atom_name,alt_notation,residue_name,residue_number,[atom coordinate]]
            x = [a[12:16].strip(),a[16],a[17:20].strip(),a[22:27].strip(),[float(a[30:38]),float(a[38:46]),float(a[46:54])]]
            if x[1] == ' ' or x[1] == 'A':
                data.append([x[0]]+x[2:])
    f.close()
    return data

def output_file_pdb(data,file_name):        #data=[atom_name,residue_name,residue+nameber,[atom coordinat]]
    tem=[]
    #print len(data)
    count=0
    for line in data:
        count+=1
        tem.append('ATOM  '+'%+5s  '%(count)+'%-4s'%(line[0])+line[1]+' A'+'%+4s'%str(line[2])+'    '+'%8.3f'%line[3][0]+'%8.3f'%line[3][1]+'%8.3f'%line[3][2]+'\n')
    file=open(file_name,'w')
    file.writelines(tem)
    file.close()

def distance(a,b):  # calculate the distance between two point
    a=np.array(a)
    b=np.array(b)
    return np.sqrt(np.dot(b-a,b-a))

def check_neib_by_group(protein,group,cluster_center, mode='vertex'):    
    '''protein=[[atom_name,res_name,res_number,[coordinate]]...]
    if use mode = 'center', must apply cluster_center. '''
    N=len(protein)
    # sort atoms by group
    atom_sorted=sorted(zip(group,protein),key=lambda x: (x[0])) 
    Group, Protein = zip(*atom_sorted)
    flag=Group[0]
    atom_name=zip(*Protein)[0]
    atom_coor=zip(*Protein)[3]
    surface=[]
    ''' center mode: choose the atom as the atom is the nearest point with the 
    center of group '''
    if mode == 'center':
        for C in cluster_center:
            group_center_atom = []
            D = 100000.0
            for atom in protein:
                d = distance(C, atom[3])
                if d <= D:
                    group_center_atom = atom
                    D = d 
            surface.append(group_center_atom)
    ''' vertex mode: choose the atom as the atom is the farthest or nearest point 
    compare with COM of protein '''
    if mode == 'vertex':
        center = np.mean(np.array(atom_coor[i]),axis=0)
        G=[]        #RAM for group information
        D=[]        #RAM for distance information
        for i in range (len(Protein)):
            if i!=len(Protein)-1:
                if Group[i]!=flag:
                    for j in range (len(G)):
                        if D[j]==max(D) or D[j]==min(D):
                            surface.append(G[j])
                            break
                    flag=Group[i]
                    G=[]
                    D=[]
                G.append(Protein[i])
                D.append(basic.distance(atom_coor[i],center))
            else:
                if Group[i]!=flag:
                    surface.append(Protein[i])
                else:
                    G.append(Protein[i])
                    D.append(basic.distance(atom_coor[i],center))
                for j in range (len(G)):
                    if D[j]==max(D) or D[j]==min(D):
                        surface.append(G[j])
                        break
    return surface

    
if __name__=='__main__':    
    ARG_VER = sys.argv
    pdbfile = sys.argv[1]
    cluster_num = int(sys.argv[2])
    protein = read_file(pdbfile)
    protein = zip(*protein)
    atom_name = protein[0]
    res_name = protein[1]
    res_number = protein[2]
    coordinate = np.array(protein[3])
    if len(coordinate)<cluster_num:
        cluster_num = len(coordinate)
    kmeans = KMeans(n_clusters=cluster_num,random_state=2).fit(coordinate)
    group = kmeans.predict(coordinate)
    centers = kmeans.cluster_centers_
    Protein = check_neib_by_group(zip(*protein),group,mode='center',cluster_center=centers)    
    if pdbfile[:3]=='pdb':
        output_file_pdb(Protein,pdbfile[3:8]+'-reduce.pdb')
    else:
        output_file_pdb(Protein,pdbfile[:5]+'-reduce.pdb')

