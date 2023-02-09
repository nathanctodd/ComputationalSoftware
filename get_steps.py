#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Creates a .xyz file for each .out file in the current directory."""

import os
import csv
from sys import argv


def main():
    """Serve as main."""
    out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.out')]
    data_fields = []
    for out_file in out_files:
        data_field = []
        with open(out_file, mode="r") as out_reader:
            count = 0
            in_vel_section = False
            for line in out_reader:
                if "  Velocities:" in line:
                    vaina = 0
                elif "###" in line or "Normal termination." in line:
                    count += 1
                elif in_vel_section and line.strip() != "":
                    vaina = 0
        data_field.append(str(out_file))
        data_field.append(count)
        data_fields.append(data_field)
	
    headers = ['file', 'steps']
    with open("num_steps.csv",'w') as file:
        #csv.write(str(DATA))
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(data_fields)
    print("")
    print("Done!")


if __name__ == "__main__":
    main()
