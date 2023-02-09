#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import readline
import numpy as np
import numpy as numpy
from re import I
from sys import argv


# -----------------  INFORMATION  -----------------
# 
#    This script is meant to parse the results from a .out file collision. It is designed with default parameters that can be changed based off of what is needed.
#    
#    INPUT:
#    It takes in a .out file from MILO software
# 
#    OUTPUT:
#    Returns in a large .csv file the speed, angle, dissociated value, and bonds broken value. 
#
#
#    RUNNING:
#    Place all of the .out files you would like to be included in the .CSV file in the same directory and then run this script using python. All of the .out files will be included.
#    It is designed around a hexane molecule hitting a wall of helium atoms, but can also take in arguments to set the parameters as necesssary. Arguments will be taken 
#    in through the command line and are outlined below.
#


print(" ------------------------------------------")
print("  _______   ___   ___   _____     ______      ")
print(" |       | |   |_|   | |     \   |   _   \    ")
print(" |       | |         | |  __  \  |  | |  |    ")
print(" |   _   | |         | | |  |  | |  |_|  |    ")
print(" |  |_|  | |         | | |__|  | |  ___  \    ")
print(" |  ___  | |  ||_||  | |      /  | |   \  \   ")
print(" |_|   |_| |__|   |__| |_____/   |_|    \__|  ")
print("                                              ")
print("                                              ")
print("    Analyical Molecular Dynamic Results       ")
print("                                              ")
print("         Created by Nathan Todd               ")
print("                                              ")
print("                                              ")

print()
print()
print()
print("Initializing...")
print()
print()
print()



# -- RESULT ORDER --
# 1. Speed
# 2. Angle
# 3. Total Steps
# 4. Dissociated Value
# 5. Bonds Broken Value

if len(argv) > 1:
    print("Arguments found.")
    if len(argv) < 9:
        print("Incorrect number of arguments: please enter the arguments in this way:")
        print("1 - # of atoms in projectile")
        print("2 - # of atoms in entire system")
        print("3 - 1st atom index in linear projectile line")
        print("4 - 2nd atom index in linear projectile line")
        print("5 - 1st atom index in  wall plane")
        print("5 - 2nd atom index in  wall plane")
        print("5 - 3rd atom index in  wall plane")
        print("6 - # of bonds that you would like to measure")
        print("7-n - atom indexes for each of the bonds. ")

    else:
        atoms_in_projectile = argv[1]
        total_atoms = argv[2]
        line_atom_1 = argv[3]
        line_atom_2 = argv[4]
        wall_atom_1 = argv[5]
        wall_atom_2 = argv[6]
        wall_atom_3 = argv[7]
        number_of_bonds = argv[8]
        bond_atoms = []
        for i in range(number_of_bonds):
            bond_atoms.append(argv[i + 9])

        print()
        print()
        print("PROJECTILE ATOMS: " + str(atoms_in_projectile))
        print("TOTAL ATOMS IN SYSTEM: " + str(total_atoms))
        print()
        print("LINE ATOM INDEX: " + str(line_atom_1))
        print("LINE ATOM INDEX: " + str(line_atom_2))
        print("WALL ATOM INDEX: " + str(wall_atom_1))
        print("WALL ATOM INDEX: " + str(wall_atom_2))
        print("WALL ATOM INDEX: " + str(wall_atom_3))
        print()
        print("NUMBER OF BONDS: " + str(number_of_bonds))
        string = ""
        for i in bond_atoms:
            string = string + str(i) + ", "
        print("BOND INDEXES: " + string)
else:
    atoms_in_projectile = 20
    total_atoms = 46


mass = np.array([12.011,1.007,1.007,1.007,12.011,1.007,1.007,\
                    12.011,1.007,1.007,12.011,1.007,1.007,\
                    12.011,1.007,1.007,12.011,1.007,1.007,1.007])/6.022e26

out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.out')]
results = []

for file in out_files:
    result = []
    result.append(file)
    with open(file, mode="r") as out_reader:
        vel_found = False
        counter = 0
        vel_matrix = []
        for line in out_reader:
            if "  Velocities " in line and vel_found == False:
                vel_found = True
            elif vel_found == True and counter < 20:
                counter += 1
                current_vel = []
                current_vel.append(float(line.strip().split()[1]))
                current_vel.append(float(line.strip().split()[2]))
                current_vel.append(float(line.strip().split()[3]))
                vel_matrix.append(current_vel)
        vel_matrix = np.array(vel_matrix)
        vel_x = np.sum((vel_matrix[:,0]*mass))/np.sum(mass)
        vel_y = np.sum((vel_matrix[:,1]*mass))/np.sum(mass)
        vel_z = np.sum((vel_matrix[:,2]*mass))/np.sum(mass)
        final_vel = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)
        final_vel = int(final_vel)
        result.append(final_vel)
    out_reader.close()
    with open(file, mode="r") as out_reader:
        xyz_found = False
        xyz_matrix = []
        counter = 0
        for line in out_reader:
            if "$molecule" in line and xyz_found == False:
                xyz_found = True
            elif xyz_found == True and len(line.strip().split()) > 2 and counter < 57:
                counter = counter + 1
                currentPos = []
                currentPos.append(float(line.strip().split()[1]))
                currentPos.append(float(line.strip().split()[2]))
                currentPos.append(float(line.strip().split()[3]))
                xyz_matrix.append(currentPos)
            if counter > 56:
                break
        wallMatrix = []
        wallMatrix.append(xyz_matrix[22])
        wallMatrix.append(xyz_matrix[26])
        wallMatrix.append(xyz_matrix[30])
        lineMatrix = []
        lineMatrix.append(xyz_matrix[4])
        lineMatrix.append(xyz_matrix[10])
        lineMatrix = np.array(lineMatrix)
        wallMatrix = np.array(wallMatrix)
        line = lineMatrix[0] - lineMatrix[1]
        QR = wallMatrix[1] - wallMatrix[0]
        QS = wallMatrix[2] - wallMatrix[0]
        normalPlane = np.cross(QR, QS)	 
        dotProduct = (normalPlane[0] * line[0]) + (normalPlane[1] * line[1]) + (normalPlane[2] * line[2])
        normalizedVainas = dotProduct / (np.linalg.norm(normalPlane) * np.linalg.norm(line))
        theta = np.arccos(normalizedVainas)
        degrees = theta * 180 / np.pi
        result.append(np.abs(degrees))
    out_reader.close()
    with open(file, mode="r") as out_reader:
        steps = -1
        for line in out_reader:
            if " Step " in line:
                steps += 1
        result.append(steps)
    out_reader.close()
    with open(file, mode="r") as out_reader:
        counter = 0
        coord_found = False
        xyz_coords = []
        xyz_temp = []
        for line in out_reader:
            if " Coordinates " in line and coord_found == False:
                coord_found = True
                xyz_temp = []
            elif coord_found == True and counter < 57:
                current_coord =[]
                counter += 1
                current_coord.append(float(line.strip().split()[1]))
                current_coord.append(float(line.strip().split()[2]))
                current_coord.append(float(line.strip().split()[3]))
                xyz_temp.append(current_coord)
            if counter > 56:
                counter = 0
                coord_found = False
                xyz_temp = np.array(xyz_temp)
                xyz_coords.append(xyz_temp)
        bond_length1 = xyz_coords[-1][0] - xyz_coords[-1][4]
        bond_length1 = np.sqrt(np.dot(bond_length1, bond_length1.T))
        bond_length2 = xyz_coords[-1][4] - xyz_coords[-1][7]
        bond_length2 = np.sqrt(np.dot(bond_length2, bond_length2.T))
        bond_length3 = xyz_coords[-1][7] - xyz_coords[-1][10]
        bond_length3 = np.sqrt(np.dot(bond_length3, bond_length3.T))
        bond_length4 = xyz_coords[-1][10] - xyz_coords[-1][13]
        bond_length4 = np.sqrt(np.dot(bond_length4, bond_length4.T))
        bond_length5 = xyz_coords[-1][13] - xyz_coords[-1][16]
        bond_length5 = np.sqrt(np.dot(bond_length5, bond_length5.T))
        resultant_bond_lengths = []
        if bond_length1 > 2.5:
            resultant_bond_lengths.append(1)
        if bond_length2 > 2.5:
            resultant_bond_lengths.append(2)
        if bond_length3 > 2.5:      
            resultant_bond_lengths.append(3)
        if bond_length4 > 2.5:      
            resultant_bond_lengths.append(4)
        if bond_length5 > 2.5:      
            resultant_bond_lengths.append(5)
        if len(resultant_bond_lengths) > 0:
            result.append(1)
        else:
            result.append(0)
        result.append(resultant_bond_lengths)
    out_reader.close()
    results.append(result)
    print("Fetching result from file: " + file)


import pandas as pd
df = pd.DataFrame(results)
df.to_csv('results.csv', index=True, header=True)
print()
print()
print()
print("Results exported to results.csv")
print("Complete!")
