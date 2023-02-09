#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Creates a .xyz file for each .out file in the current directory."""

import os
from re import I
from sys import argv


def main():
    """Serve as main."""
    out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.out')]

    for out_file in out_files:
        step_number = 0
        with open(out_file, mode="r") as out_reader:
            num_atoms = None
            final_vel_lines = list()
            current_vel = list()
            if (len(argv) > 1):
                num_atoms = argv[1]
            count = 0
            atom_number = 0
            in_vel_section = False
            for line in out_reader:
                if "Step" in line:
                    # print("step found")
                    step_number += 1
                if "  Velocities:" in line:
                    in_vel_section = True
                    atom_number = 0
                elif "###" in line or "Normal termination." in line:
                    if count > 5:
                        if num_atoms is None:
                            num_atoms = str(len(current_vel))
                        if len(current_vel) > 0:
                            final_vel_lines.append(num_atoms)
                            final_vel_lines.append("")
                            final_vel_lines.extend(current_vel)
                        current_vel = list()
                        in_vel_section = False
                    count += 1
                elif in_vel_section and line.strip() != "":
                    if (atom_number < int(num_atoms)):
                        current_vel.append(line.strip())
                    atom_number += 1


        with open(out_file[:-4] + ".vel", mode="w") as vel_writer:
            vel_writer.write(str(step_number - 1) + "\n")
            for line in final_vel_lines:
                vel_writer.write(line + "\n")
        print(out_file[:-4] + ".vel")
    print("Done!")
        


if __name__ == "__main__":
    main()
