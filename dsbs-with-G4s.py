# Author: Ella Chee
# Date: August 2024
# Description: Script to create .bed files of DSBs with G4s, 
# by intersecting DSBs with corresponding G4Catchall data.

# Imports
from pybedtools import BedTool 
import pybedtools
import pandas as pd
import os
import csv

# Convert name to G4catchall file name
def g4catchall_func(file):
    '''Converts DSB file name to G4catchall file name'''
    return file.replace('_remove_neg_starts', 'converted_valid')

# Convert name to shuffle file name
def shuffle_func(file):
    '''Converts DSB file name to shuffle file name'''
    return file.replace('.breakends', '.shuffle.breakends')

# File paths and names
all_dsb_path = '/work/daylab/fall-2024-coops/INDUCE-seq-beds/all_breaks_clean/'
g4catchall_path = '/work/daylab/fall-2024-coops/INDUCE-seq-beds/g4_catchall_clean/'
g4_dsbs_path = '/work/daylab/fall-2024-coops/INDUCE-seq-beds/g4_dsbs/'

dsb_files = ['siBAZ2A_BRACO19_18h.breakends_ext100_remove_neg_starts.bed',
             'siBAZ2A_veh_0h.breakends_ext100_remove_neg_starts.bed',
             'siBAZ2B_BRACO19_18h.breakends_ext100_remove_neg_starts.bed', 
             'siBAZ2B_veh_0h.breakends_ext100_remove_neg_starts.bed',
             'siControl_BRACO19_18h.breakends_ext100_remove_neg_starts.bed', 
             'siControl_veh_0h.breakends_ext10_remove_neg_starts.bed']

# Create list of G4catchall file names
g4catchall_data = [g4catchall_func(file) for file in dsb_files]

# Intersection to remove DSBs for single file in list
def remove_dsbs(i):
    '''Intersect all DSBs with G4Catchall data to create dataset of DSBs without G4s'''
    d = dsb_files[i]
    g = g4catchall_data[i]
    dsb_file = all_dsb_path + d
    g4catchall_file = g4catchall_path + g

    dsb_bed = BedTool(dsb_file).sort()
    g4catchall_bed = BedTool(g4catchall_file).sort()

    # Intersection of DSBs and G4catchall to get DSBs w/ G4s
    intersection = dsb_bed.intersect(g4catchall_bed, v=True)
    output_file = g4_dsbs_path + d.replace('.bed', '_G4_DSBs.bed')
    intersection.moveto(output_file)
    print(f"Saved intersection to {output_file}")

# Remove DSBs for all files in list
for i in range(len(dsb_files)):
    remove_dsbs(i)
print('All DSBs with G4s saved to g4_dsbs folder')