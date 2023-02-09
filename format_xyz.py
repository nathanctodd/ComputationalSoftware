
import readline
import numpy as np
from sys import argv
import os



allPercentages = []



xyz_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.xyz')]
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
               step.append("\n\n")
               num_steps = num_steps + 1
           else:
               if len(line) > 2:
                   step.append(line)

    with open((fileName.split(".")[0] + "_formated.xyz"), "a") as out:
        for step in steps:
            for line in step:
                out.write(line)
            #out.write("\n")
                

print("Done!")
# with open("kinetic_energies.csv",'w') as csv:
# 	csv.write(allPercentages)
# print("_distances.csv file written")

#import pandas as pd
#df = pd.DataFrame(allPercentages)
#df = df.T
#df.to_csv('_distances.csv', index=True, header=True)
