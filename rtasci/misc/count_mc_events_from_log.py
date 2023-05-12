# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Leonardo Baroncelli <leonardo.barconcelli@inaf.it>
# *******************************************************************************

import os 
import argparse

if __name__=='__main__':

    # get input
    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--input-dir', type=str, help='input dir')
    args = parser.parse_args()

    background_counts = 0
    events_counts = 0

    for root, dirs, files in os.walk(args.input_dir):
        for f in files:
            if ".log" in f and ".skymap" not in f:
                f = os.path.join(root, f)
                # print(f)
                with open(f, "r") as fh:
                    lines = fh.readlines()  
                    for line in lines:
                        if ":  MC source events ..........: " in line and "[run0406_ID000126]" in line:
                            #print("\n->",line)
                            line = line.split("MC source events ..........: ")[1].strip().replace("\n","")   
                            #print("line:",line)
                            counts = line.split(" [run0406_ID000126]")[0]
                            #print("counts: ", counts)
                            counts = int(counts)
                            #print(f"source counts: {counts}")
                            events_counts += counts
                        if ":  MC background events ......: " in line:
                            #print("\n->",line)
                            line = line.split("MC background events ......: ")[1].strip().replace("\n","")
                            #print(f"line:{line}")
                            counts = int(line)
                            #print("counts: ", counts)
                            #print(f"bkg counts: {counts}")
                            background_counts += counts
    
    print("Background MC simulated events:",background_counts)
    print("Source MC simulated counts:",events_counts)
    print("Total MC simulated events:",background_counts+events_counts)
