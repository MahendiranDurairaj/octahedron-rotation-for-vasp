# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:43:56 2019

@author: jliubt
"""

# rotate the LiCoO2 primitive cell to align the Co-O bond in octahedra with inner coordiante
# to accurate split the orbital
# 
# IMPORTANT NOTES TO RUN THIS CODE!!!
# 1) The Co atom is shifted to [0, 0, 0]
# 2) Therefore, the unit axis in rotation matrix can be determined by the Co-O bond
# 3) Only cartesian coordinate is considered 
import numpy as np
from numpy.linalg import norm
import sys

# define the rotation matrix according to two different coordinate
# https://www.continuummechanics.org/rotationmatrix.html
def Rotation_matrix(coordinate_prime, coordinate_new):
    
    R = np.zeros_like(coordinate_prime)
    for i in range(3):
        for j in range(3):
            R[i,j] = np.dot(coordinate_prime[j], coordinate_new[i])
            
    return R

# read in the POSCAR
def read_POSCAR(filename):
    
    with open(filename, 'r') as f:
        sys = f.readline().rstrip('\n')
        scale = float(f.readline().rstrip('\n'))
        lattice = np.zeros((3,3))
        for i in range(3):
            line = f.readline().rstrip('\n').split()
            lattice[i] = np.array([float(s)*scale for s in line])
        # write atom type and atom number
        atom = {}
        atom["element"] = f.readline().rstrip('\n').split()
        atom["number"] = np.array([int(s) for s in f.readline().rstrip('\n').split()])
        atom_number = sum(atom["number"])
        # skip the line
        f.readline()
        # IMPORTANT INFORMATION!!!
        # position in cartesian coordinate
        pos = np.zeros((atom_number,3))
        for i in range(atom_number):
            line = f.readline().rstrip('\n').split()
            pos[i] = np.array([float(s) for s in line])
        
        return sys, lattice, atom, pos

# write the POSCAR
def write_POSCAR(sys, atom, lattice, pos, filename):
    
    atom_type = len(atom["element"])
    atom_number = sum(atom["number"])
    
    with open(filename, 'w') as f:
        f.write("{}  \n1.0\n".format(sys))
        for i in range(3):
            f.write("    {0:.10f}    {1:.10f}    {2:.10f}\n".format(lattice[i,0], lattice[i,1], lattice[i,2]))
       # write atom type and atom number
        for i in range(atom_type):
            f.write("  {}".format(atom["element"][i]))
        f.write("\n")
        for i in range(atom_type):
            f.write("  {}".format(atom["number"][i]))
        f.write("\nCartesian\n")
        # IMPORTANT INFORMATION!!!
        # write position in cartesian coordinate
        for i in range(atom_number):
            f.write("    {0:.9f}    {1:.9f}    {2:.9f}\n".format(pos[i,0], pos[i,1], pos[i,2]))

# read in the .xyz file
def read_xyz(filename):
    
    with open(filename, 'r') as f:
        atom_number = int(f.readline().rstrip('\n'))
        sys = f.readline().rstrip('\n')
        element = []
        pos = np.zeros((atom_number,3))
        
        for i in range(atom_number):
            line = f.readline().rstrip('\n').split()
            element.append(line[0])
            pos[i] = np.asarray([float(line[i+1]) for i in range(3)])
            
        return atom_number, sys, element, pos
     
# write the .xyz file
def write_xyz(atom_number, element, sys, pos, filename):
    
    with open(filename, 'w') as f:
        f.write("{}\n{}\n".format(atom_number, sys))
        
        for i in range(atom_number):
            f.write("{0}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(element[i], pos[i,0], pos[i,1], pos[i,2]))


if __name__=="__main__":
    filename = sys.argv[1]
    # original coordiante
    coordinate_prime = np.asarray([[1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    # IMPORTANT PART!!!
    # new coordinate, should define for different system
    z_new = np.asarray([[0.710588, 0.209747, 1.783932]])
    z_new = z_new/norm(z_new)
    x_new = np.asarray([[1.510614, -1.027360, -0.627628]])
    y_new = np.cross(z_new, x_new)/norm(np.cross(z_new, x_new))
    x_new = np.cross(y_new, z_new)
    coordinate_new = np.concatenate((x_new, y_new, z_new))
    
    # calculate the rotation matrix
    R = Rotation_matrix(coordinate_prime, coordinate_new)
    
    # filename = "POSCAR"
    # read in the POSCAR
    _, lattice, atom, pos = read_POSCAR("{}.vasp".format(filename))
    lattice_new = np.matmul(lattice, R.T)
    pos_new = np.matmul(pos, R.T)
    # write the POSCAR
    sys = "rotated cell"
    write_POSCAR(sys, atom, lattice_new, pos_new, "{}_rotated.vasp".format(filename))
    
    # read in the .xyz file
    atom_number, _, element, pos = read_xyz("{}.xyz".format(filename))
    pos_new = np.matmul(pos, R.T)
    # write the .xyz file
    write_xyz(atom_number, element, sys, pos_new, "{}_rotated.xyz".format(filename))
