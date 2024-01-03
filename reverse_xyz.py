
import readline
import numpy as np
from sys import argv
import os


class Reverser():
    
    def reverse_xyz(xyz_files):
        for fileName in xyz_files:
            steps = []
            num_steps = 0
            numOfAtoms = 0
            print(fileName)
            with open(fileName,'r') as f:
                file_iter = iter(f)
                step=[]
                for line in file_iter:
                    if len(line.strip().split()) == 1:
                        numOfAtoms = line.strip().split()[0]
                        if len(step) > 0:
                            steps.append(step)
                        step = []
                        step.append(str(numOfAtoms))
                        num_steps = num_steps + 1
                    else:
                        step.append(line)

            with open((fileName.strip(".") + "_reversed.xyz"), "a") as out:
                for step in reversed(steps):
                    out.write("\n")
                    for line in step:
                        out.write(line)

xyz_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.xyz')]
reverser = Reverser()
reverser.reverse_xyz(xyz_files)
                
print("Done!")

