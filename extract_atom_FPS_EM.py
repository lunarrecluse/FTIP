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
from scipy.spatial import distance_matrix
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

def extract_point_label(distance_matrix, label, point_num, mode = 'MAX'):    # extract point based on geometry distances
    dist_M = np.array(distance_matrix)
    if len(dist_M) <= point_num:
        point_num = len(dist_M)
    EP = []
    EP_label = []
    if mode == 'MAX':
        # step 1: pick initial 2 points
        x = np.where(dist_M == np.max(dist_M))
        EP.append(x[0][0])
        EP.append(x[1][0])
        EP_label.append(label[x[0][0]])
        EP_label.append(label[x[1][0]])
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
            EP_label.append(label[x])
            sum_list = sum_list + dist_M[x]
    if mode == 'MIN':
        gamma = np.mean(dist_M) / 100000.0
        dist_M = 1 / (dist_M + gamma)
        # step 1: pick initial 2 points
        x = np.where(dist_M == np.min(dist_M))
        EP.append(x[0][0])
        EP.append(x[1][0])
        EP_label.append(label[x[0][0]])
        EP_label.append(label[x[1][0]])
        sum_list = dist_M[x[0][0]] + dist_M[x[1][0]]
        # step 2: pick rest points until reach the desired amount
        while len(EP) < point_num:
            x = np.where(sum_list == np.min(sum_list))[0][0]
            if x in EP:
                print "ERROR, point exist!"
            EP.append(x)
            EP_label.append(label[x])
            sum_list = sum_list + dist_M[x]
    return EP, EP_label

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

def get_points_byLabel(coor, label, target_label):
    result = []
    result_label = []
    for i in range (len(label)):
        if label[i] in target_label:
            result.append(coordinate[i])
            result_label.append(label[i])
    return result, result_label


def Label_point(coor,window):
    x,y,z = np.int_(np.array(coor)/window)+1
    return str(x)+'|'+str(y)+'|'+str(z)

def label_by_window(coor,window):  
    coor = shift(coor)
    center = {}     # {label: [[coor],[count]]}
    label = ['000' for i in range (len(coor))]
    for i in range (len(coor)):
        label[i] = Label_point(coor[i],window)
        if label[i] in center.keys():
            center.update({label[i]:[center[label[i]][0]+coor[i], center[label[i]][1]+1]})
        else:
            center.update({label[i]:[coor[i],1]})
    return label, center

def shift(coor,origin=[0.0,0.0,0.0]):
    coor = np.array(coor)
    return coor - np.min(coor,axis=0)
    
if __name__=='__main__':    
    ARG_VER = sys.argv
    pdbfile = sys.argv[1]
    point_num = int(sys.argv[2])
    grid_size = int(sys.argv[3])
    Mode = sys.argv[4]
    read_begin = 0
    coordinate = []
    with open (pdbfile) as f:
        for line in f:
            if line == "      <Coordinate point='\n":
                read_begin = 1 
                continue
            if line == "      '/>\n" and read_begin == 1:
                break
            if read_begin == 1:
                coordinate.append(map(float,line.split()))
    label, center = label_by_window(coordinate, grid_size)
    Center = {}
    for item in center:
        Center.update({item:center[item][0]/center[item][1]})
    center_labels, center_points = zip(*[[k,v] for k,v in Center.items()])
    # first step extract labels
    DM = distance_matrix(center_points, center_points)
    points_1_index,picked_label = extract_point_label(DM, center_labels, point_num, mode = Mode)
    # second step: get poins with such labels
    points_2, label_2 = get_points_byLabel(coordinate, label, picked_label)
    # third step: final extract
    DM = distance_matrix(points_2,points_2)
    points_final_index, points_final_label = extract_point_label(DM, label_2, point_num, mode = Mode)
    points_final = [points_2[x] for x in points_final_index]
    data = [['X', 'YYY', 'Z', x ] for x in points_final]
    output_file_pdb(data,pdbfile[0:5]+'-'+Mode+'.pdb')        #data=[atom_name,residue_name,residue+nameber,[atom coordinat]]

