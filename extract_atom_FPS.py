#!/usr/bin/env python
#(c) Yuan Zhang in Department of Statistic, FSU
# Nov.1, 2017
# This script is used to get the key atoms on protein surface

import sys
import numpy as np
import pylab as py
import matplotlib
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift
import copy

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

def distance(a,b):  # calculate the distance between two point
    a=np.array(a)
    b=np.array(b)
    return np.sqrt(np.dot(b-a,b-a))

def dist_matrix(data):  # calculate the distance matrix of data, symmetric matrix with diagonal=0
    n = len(data)
    matrix = np.array([[0.0 for i in range (n)] for j in range (n)])
    for i in range (n-1):
        for j in range (i+1,n): 
            D = distance(data[i],data[j])
            matrix[i][j] = D 
            matrix[j][i] = D 
    return matrix

def extract_point(distance_matrix, point_num, mode = 'MAX'):    # extract point based on geometry distances
    dist_M = np.array(distance_matrix)
    if len(dist_M) <= point_num:
        point_num = len(dist_M)
    EP = []
    if mode == 'MAX':
        # step 1: pick initial 2 points
        x = np.where(dist_M == np.max(dist_M))
        EP.append(x[0][0])
        EP.append(x[1][0])  
        sum_list = dist_M[x[0][0]] + dist_M[x[1][0]]
        # step 2: pick rest points until reach the desired amount
        while len(EP) < point_num:
            tem_sum_list = copy.deepcopy(sum_list)
            for item in EP: 
                tem_sum_list[item] = 0 
            x = np.where(sum_list == np.max(tem_sum_list))[0][0]
            if x in EP: 
                print "ERROR, point exist!"
                print np.where(sum_list == np.max(Sum_list))
            EP.append(x)
            sum_list = sum_list + dist_M[x]
    if mode == 'MIN':
        gamma = np.mean(dist_M) / 100000.0
        dist_M = 1 / (dist_M + gamma)
        # step 1: pick initial 2 points
        x = np.where(dist_M == np.min(dist_M))
        EP.append(x[0][0])    
        EP.append(x[1][0])
        sum_list = dist_M[x[0][0]] + dist_M[x[1][0]]
        # step 2: pick rest points until reach the desired amount       
        while len(EP) < point_num:
            x = np.where(sum_list == np.min(sum_list))[0][0]
            if x in EP: 
                print "ERROR, point exist!"
            EP.append(x)
            sum_list = sum_list + dist_M[x] 
    return EP    

def output_file_pdb(data,file_name):        #data=[atom_name,residue_name,residue+nameber,[atom coordinat]]
    tem=[]
    print len(data)
    count=0
    for line in data:
        count+=1
        tem.append('ATOM  '+'%+5s  '%(count)+'%-4s'%(line[0])+line[1]+' A'+'%+4s'%str(line[2])+'    '+'%8.3f'%line[3][0]+'%8.3f'%line[3][1]+'%8.3f'%line[3][2]+'\n')
    file=open(file_name,'w')
    file.writelines(tem)
    file.close()


	
if __name__=='__main__':	
	ARG_VER = sys.argv
	pdbfile = sys.argv[1]
	point_num = int(sys.argv[2])
	mode = sys.argv[3]
	protein = read_file(pdbfile)
	protein = zip(*protein)
	atom_name = protein[0]
	res_name = protein[1]
	res_number = protein[2]
	coordinate = np.array(protein[3])
	distance_matrix = dist_matrix(coordinate)
	EP = extract_point(distance_matrix, point_num, mode = mode)
	Protein = []
	for i in EP:
		Protein.append([atom_name[i], res_name[i], res_number[i], coordinate[i]])
	output_file_pdb(Protein,pdbfile[0:5]+'-'+mode+str(point_num)+'.pdb')
