import re
import numpy as np
import string
from pylammpsmpi import LammpsLibrary
import MDAnalysis as mda
import MDAnalysis.analysis.align
from MDAnalysis import *
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF,DCD

def load_files(fname, debug=False):
    load_file_export = open(fname, 'r')
    lines = load_file_export.readlines()
    no_lines = len(lines)
    frames_cnt = 0
    format_data= None
    no_bonds = 0
    bonds = 0
    for i in range(no_lines):
        line = lines[i]
        if 'item: number of atoms'.upper() in line:
            format_data='LAMMPSDUMP'
            break
        i=i+1

    frames_cnt = 0
    load_file = MDAnalysis.Universe(fname,format=format_data)
    if format_data == None:
        bonds = load_file.bonds.to_indices()
        no_bonds = len(load_file.bonds)
    return load_file,no_bonds
