#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import readline
import numpy as np
import numpy as numpy
from re import I
from sys import argv


#     -------------------------------------------------
#       _______   _______    ___________    _______       
#      |       | |   __   \ |           |  /        \     
#      |       | |  |  |  | |___     ___| |    __    |    
#      |   _   | |  |__|  |     |   |     |   |  |   |    
#      |  |_|  | |   ___  \     |   |     |   |__|   |    
#      |  ___  | |  |   \  \    |   |     |          |    
#      |_|   |_| |__|    \__|   |___|      \________/     
#                                                  
#                                                 
#        Analytical Resultant Trajectory Operation        
#                                                     
#                Created by Nathan Todd                      
#                                                   
#
#
# -----------------  INFORMATION  -----------------
# 
#    This script is meant to parse the results from a .out file collision. It is designed with default parameters that can be changed based off of what is needed.
#    
#    INPUT:
#    There are three different working pieces of software combined into ARTO. 
#
#    [Job Index] 
#       0 - Hexane Sampling
#       1 - Trajectory Build & Submit
#       2 - Result Parser
#
#   1. Hexane Sampling - called by using the following command:
#
#       python  Analytical_Resultant_Trajectory_Operation.py  [Job Index]  [# of Hexane to Sample]
#   
#       - Hexane Sampling takes the "hexane.in" file and submits it to Milo to run however many times is needed. 
#       - This is to get enough vibrational sampling data for the trajectories that are run in the next step.
#       - Hexane Sampling relies on a "hexane.in" and "hexane.sh" to be able to generate and submit all of the necessary hexane sampling runs.
#       - Note: wait for all hexane sampling files to finish before proceeding to step 2.
#
# 
#
#   2. Trajectory Build & Submit - called by using the following command:
#       
#       python Analytical_Resultant_Trajectory_Operation.py  [Job Index] [Velocity to add in x-direction (m/s)]
#
#       - This takes all of the previously generated hexane sampling .out files, locates the velocities from them, then adjusts them according to the needed velocity to add.
#       - After this, the code then puts them into a newly copied .in file and generates a .sh file to help submit them. It does this automatically too. 
#       - The Trajectory Build & Submit relies on a "hexane_trajectory_.in" and a "hexane_trajectory_.sh" to repeatedly copy, insert, and run each trajectory.
#       - Note, this part of the code also cleans up some of the slurm files, as well as organizes the hexane sampling. Also please wait until trajectories are all finished before proceeding to step 3.
# 
# 
#
#   3. Result Parser - called by using the following command:
#       
#       python Analytical_Resultant_Trajectory_Operation.py  [Job Index] 
#
#       - This part of the code takes all of the new trajectory .out files and creates a .csv file full of all data. This data includes:
#
#           # -- RESULT ORDER --
#               1. Speed
#               2. Angle
#               3. Total Steps
#               4. Dissociated Value
#               5. Bonds Broken Value
#
#       - After this, the code then parses all of the .out files and generates the kinetic energy graph, returned as kinetic_energies.csv.
#       - Note: This produces two .csv files, both of which include rows/columns for every single trajectory. This allows for easy comparison and graphing, as well as consolidation.
# 
#
#
#    SYSTEM REQUIREMENTS:
#       - In order for ARTO to work properly, this python script must have access to a supercomputer using SLURM, as well as MILO created by the Ess Group at Brigham Young University
#       - Also, the proper template files must be included for this script to work properly at the end. 
#
#
#
#

class Analytical_Resultant_Trajectory_Operation():

    def write_velocity_files(out_files):
        total_number = 0
        for out_file in out_files:
            total_number += 1
            with open(out_file, mode="r") as out_reader:
                velocities = []
                found_sampling = False
                for line in out_reader:
                    if "  Energy Summary (kcal/mol):" in line:
                        found_sampling = False
                    elif found_sampling == True:
                        velocities.append(line)
                    elif "### Velocity Initialization Summary" in line:
                        found_sampling = True
            vel_file_name = "vel" + str(total_number) + ".velocity"
            velocities.pop(0)
            with open(vel_file_name, mode="w") as out_writer:
                for line in velocities:
                    out_writer.write(line)
        return total_number
                    
    def write_vel_files(velocity_files):
        for velocity_file in velocity_files:
            with open(velocity_file, mode="r") as velocity_reader:
                final_lines = []
                for line in velocity_reader:
                    split_line = line.split()
                    vel1 = float(split_line[1]) + int(velocity_to_add)
                    vel2 = float(split_line[2])
                    vel3 = float(split_line[3])
                    new_line = split_line[0] + "  " + str(vel1) + "  " + str(vel2) + "  " + str(vel3) + "  \n"
                    final_lines.append(new_line)
            with open(velocity_file.split('.')[0] + ".vel", mode="w") as velocity_writer:
                for line_given in final_lines:
                    velocity_writer.write(line_given)

class Trajectory_Operator():
    
    def write(out_files):
        for out_file in out_files:
            step_number = 0
            with open(out_file, mode="r") as out_reader:
                num_atoms = None
                final_vel_lines = list()
                current_vel = list()
                if (len(argv) > 1):
                    num_atoms = 57
                count = 0
                atom_number = 0
                in_vel_section = False
                for line in out_reader:
                    if "Step" in line:
                        # print("step found")
                        step_number += 1
                    if "  Velocities (meter/second):" in line:
                        in_vel_section = True
                        atom_number = 0
                    elif "Energy Summary" in line or "Normal termination." in line:
                        # if count > 5:
                            if num_atoms is None:
                                num_atoms = str(len(current_vel))
                            if len(current_vel) > 0:
                                final_vel_lines.append(num_atoms)
                                final_vel_lines.append("\n")
                                final_vel_lines.extend(current_vel)
                            current_vel = list()
                            in_vel_section = False
                        # count += 1
                    elif in_vel_section and line.strip() != "":
                        if (atom_number < int(num_atoms)):
                            current_vel.append(line)
                        atom_number += 1


            with open(out_file[:-4] + ".vel_file", mode="w") as vel_writer:
                vel_writer.write(str(step_number - 1) + "\n")
                for line in final_vel_lines:
                    vel_writer.write(str(line))
            print(out_file[:-4] + ".vel_file")
            
    def vel_files(vel_files):
        allPercentages = []
        for fileName in vel_files:
            print("Currently working on file: " + fileName)
            dataSteps = []
            num_steps = 0
            numOfAtoms = 0
            with open(fileName,'r') as f:
                # data = f.read()
                num_steps = f.readline().strip()
                num_steps = int(num_steps, 10)
                for i in range(num_steps):
                    numOfAtoms = int(f.readline().strip(), 10)
                    #space = f.readline()
                    #print(numOfAtoms)
                    data_step = []
                    for atomIndex in range(numOfAtoms):
                        line = f.readline()
                        data_step.append(line)
                    dataSteps.append(data_step)



            new_dataSteps = []
            for i in dataSteps:
                temp = []
                for j in range(20):
                    temp.append(i[j])
                new_dataSteps.append(temp)


            total_kinetic_energies = []
            total_translational_energies = []
            percentages = []
            percentages2 = []
            givenVelocityPercentages = []


            total_kinetics = []

            for dataStep in new_dataSteps:
                velocities = np.zeros((20,3))
                for i in range(numOfAtoms):
                    if i < 20:
                        velocities[i, 0] = float(dataStep[i].split()[1])
                        velocities[i, 1] = float(dataStep[i].split()[2])
                        velocities[i, 2] = float(dataStep[i].split()[3])
            
                #Setup mass in KG
                mass = np.array([12.011,1.007,1.007,1.007,12.011,1.007,1.007,\
                        12.011,1.007,1.007,12.011,1.007,1.007,\
                        12.011,1.007,1.007,12.011,1.007,1.007,1.007])/6.022e26

                # Calculate total kinetic energy for dataStep
                KE_0 = np.sum(mass/2*np.sum(velocities**2,axis=1))
                #print(KE_0)
                total_kinetic_energies.append(KE_0)

                # Calculate the center of mass velocities
                vel_x = np.sum((velocities[:,0]*mass))/np.sum(mass)
                vel_y = np.sum((velocities[:,1]*mass))/np.sum(mass)
                vel_z = np.sum((velocities[:,2]*mass))/np.sum(mass)

                # Calculate the velocities with the COM velocity subtracted out
                vels = np.zeros((20,3))
                vels[:,0] = velocities[:,0]-vel_x
                vels[:,1] = velocities[:,1]-vel_y
                vels[:,2] = velocities[:,2]-vel_z

                # The kinetic energy minus the translational component (So just vibrational and rotational)
                # Kinetic energy without translational component
                total_kinetic_energy_minus_translational = np.sum(mass/2*np.sum(vels**2,axis=1))
                percent3 = round(np.abs(total_kinetic_energy_minus_translational)/KE_0*100,1)
                total_kinetics.append(percent3)

                # Print the % of kinetic energy not in vibrations

            allPercentages.append(total_kinetics)
        return allPercentages


if __name__ == "__main__":
        

    print(" -------------------------------------------------")
    print("  _______   _______    ___________    _______       ")
    print(" |       | |   __   \ |           |  /        \     ")
    print(" |       | |  |  |  | |___     ___| |    __    |    ")
    print(" |   _   | |  |__|  |     |   |     |   |  |   |    ")
    print(" |  |_|  | |   ___  \     |   |     |   |__|   |    ")
    print(" |  ___  | |  |   \  \    |   |     |          |    ")
    print(" |_|   |_| |__|    \__|   |___|      \________/     ")
    print("                                                     ")
    print("                                                     ")
    print("   Analytical Resultant Trajectory Operation        ")
    print("                                                     ")
    print("         Created by Nathan Todd                      ")
    print("                                                     ")
    print("                                                     ")


    print()
    print()
    print()
    print("Initializing...")
    print()
    print()
    print()
    print()
    print()




    #
    #
    # 1. Sampling Hexane
    #
    #

    stage = int(argv[1])
    if stage == 0:
        iterations = argv[2]
        print("Running " + str(iterations) + " iterations...")
        print()
        print()
        print()
        os.system("for i in {1.." + str(iterations) + "}; do sbatch hexane_sampling.sh; done")
        print()
        print()
        print()
        print("Jobs have been submitted for hexane sampling:")
        print("Please wait for sampling to complete on super computer before proceeding.")
        print()
        print()
        print()



#
#
# 2. Build & Run Trajectories
#
#

    elif stage == 1:
        
        ARTO = Analytical_Resultant_Trajectory_Operation()
        
        os.system("rm *.xyz")
        #Add parameter to include random velocity, ranged velocity, or set velocity
        velocity_to_add = argv[2]
        out_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane') and f.endswith('.out')]
        print("Fetching results: " + str(len(out_files)) + " sampling files found. Creating trajectories:\n\n")
        total_number = ARTO.write_velocity_files(out_files)
        print("Velocities found and recorded\n\n")
        
        
        velocity_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.velocity')]
        ARTO.write_vel_files(velocity_files)
        print("Velocities added to each .vel file")
        print("Deleting old velocities...")
        os.system("rm *.velocity")
        print('\n\n')
        
        
        submit_file_lines_read = []
        trajectory_in_file_lines_read = []
        readable_submit_file = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane_trajectory_.sh')]
        with open(readable_submit_file[0], mode="r") as submit_readable:
            for line in submit_readable:
                submit_file_lines_read.append(line)
        trajectory_files_ = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane_trajectory_.in')]
        with open(trajectory_files_[0], mode="r") as trajectory_readable:
            for line in trajectory_readable:
                trajectory_in_file_lines_read.append(line)

        print("Creating .in and .sh files for submission...")
        print()
        for i in range(1, total_number + 1):
            os.system("cp hexane_trajectory_.in hexane_trajectory_" + str(i) + ".in")
            os.system("cp hexane_trajectory_.sh hexane_trajectory_" + str(i) + ".sh")

        trajectory_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane_trajectory_') and f.endswith('in')]
        trajectory_submit_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane_trajectory_') and f.endswith('sh')]

        for i in range(1, total_number + 1):
            trajectory_submit_file = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith ('hexane_trajectory_' + str(i) + ".sh")]
            with open(trajectory_submit_file[0], mode="w") as trajectory_submit_writer:
                for index in submit_file_lines_read:
                    if "export JOB_NAME=" in index:
                        trajectory_submit_writer.write("export JOB_NAME=hexane_trajectory_" + str(i) + "\n")
                    else:
                        trajectory_submit_writer.write(index)
        vel_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith(".vel")]
        for i in range(1, total_number + 1):
            trajectory_in_file = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith ('hexane_trajectory_' + str(i) + ".in")]
            vel_lines = []
            with open(vel_files[i - 1], mode="r") as vel_reader:
                for line in vel_reader:
                    vel_lines.append(line)
            with open(trajectory_in_file[0], mode="w") as trajectory_in_writer:
                found_vel_section = False
                current_index = 0
                for line in trajectory_in_file_lines_read:
                    if "$velocities" in line:
                        found_vel_section = True
                        trajectory_in_writer.write(line)
                        current_index = 0
                    elif current_index >= 20 and found_vel_section == True:
                        found_vel_section = False
                        trajectory_in_writer.write(line)
                    elif found_vel_section == True:
                        trajectory_in_writer.write(vel_lines[current_index].split()[1] + " " + vel_lines[current_index].split()[2] + " " + vel_lines[current_index].split()[3] + "\n")
                        current_index += 1
                    elif found_vel_section == False:
                        trajectory_in_writer.write(line)

        print(".in and .sh files are completed and ready for submission.")
        print("Submitting to super computer:")
        submitable_sh_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('hexane_trajectory_') and f.endswith('.sh')]
        for i in submitable_sh_files:
            print("Submitting " + str(i))
            os.system("sbatch " + i)
        print()
        print("Done! Please wait for all trajectories to finish before using the result parser.")






    #
    #
    # 3. Result Parser
    #
    #

    else:
        print("Running result parser:")
        print()
        print()
        print()
        print()
        print("Storing slurm and hexane sampling files:")
        os.system("mkdir slurm_files")
        os.system("mv slurm-* slurm_files")
        os.system("mkdir hexane_sampling")
        os.system("mv hexane_5* hexane_sampling")
        atoms_in_projectile = 20
        total_atoms = 57
        print()
        print("Done")
        print()
        print()


        mass = np.array([12.011,1.007,1.007,1.007,12.011,1.007,1.007,\
                            12.011,1.007,1.007,12.011,1.007,1.007,\
                            12.011,1.007,1.007,12.011,1.007,1.007,1.007])/6.022e26

        out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                        f.endswith('.out') and f.startswith('hexane_trajectory')]
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
                    elif xyz_found == True and len(line.strip().split()) > 2 and counter < total_atoms:
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
                    elif coord_found == True and counter < total_atoms:
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

        Traj_Operator = Trajectory_Operator()
        
        Traj_Operator.write(out_files)

        
        print("Done parsing initial results!")
        print()
        print()
        print("Continuing on to parse velocities and get kinetic energy information:")
        print()
        print()

        vel_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                    f.endswith('.vel_file')]
        
        allPercentages = Traj_Operator.vel_files(vel_files)
        

        import pandas as pd
        df = pd.DataFrame(allPercentages)
        df = df.T
        df.to_csv('kinetic_energies.csv', index=True, header=True)
        print("Done!")



        import pandas as pd
        df = pd.DataFrame(results)
        df.to_csv('results.csv', index=True, header=True)
        print()
        print()
        print()
        print("Results exported to results.csv")
        print("Complete!")


    print()
    print()
    print("Thank you for using ARTO!")
    print()
    print()
