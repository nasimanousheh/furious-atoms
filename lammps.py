import re
import numpy as np
import string


def load_lammps(fname, debug=False):

    lammps_file = open(fname, 'r')
    lines = lammps_file.readlines()
    no_lines = len(lines)
    lammps_dix = {'index': {}}
    frames_cnt = 0

    for i in range(no_lines):
        line = lines[i]
        if 'item: number of atoms'.upper() in line:
            no_atoms = int(lines[i+1])
            lammps_dix['index'][frames_cnt] = {'box': [],
                                               'coords': [],
                                               'no_atoms': no_atoms, 'bonds': []}
        if 'box'.upper() in line:
            # line_lx = list(map(float, lines[i+1].split('\n')[0].split('\t')))
            # line_ly = list(map(float, lines[i+2].split('\n')[0].split('\t')))
            # line_lz = list(map(float, lines[i+3].split('\n')[0].split('\t')))
            box_lx = 50 #(np.abs(line_lx[0])+np.abs(line_lx[1]))
            box_ly = 50 #(np.abs(line_ly[0])+np.abs(line_ly[1]))
            box_lz = 50 #(np.abs(line_lz[0])+np.abs(line_lz[1]))
            # print( box_lx, box_ly, box_lz)
            lammps_dix['index'][frames_cnt]['box'] = [box_lx, box_ly, box_lz]

        if 'item: atoms'.upper() in line:
            words = lines[i].split()
            try:
                item_index = words.index('index')
                x_index = words.index('x')
                y_index = words.index('y')
                z_index = words.index('z')
                atom_type_index = words.index('type')
                no_columns = 5
                if debug:
                    print('Standard x y z seen')
            except ValueError:
                no_columns = len(lines[i + 1].split())
                if debug:
                    print('Non standard x y z seen')
            frame_as_list_text = lines[i+1: i+1 + no_atoms]
            frame_positions = np.genfromtxt(frame_as_list_text)
            lammps_dix['index'][frames_cnt]['coords'] = frame_positions
            # frames_cnt += 1

        no_bonds = 4880
        if 'Bonds' in line:
            print("bond is found")
            frame_as_list_text_for_bonds = lines[i+1: i+1 + no_bonds]
            bond_locations = np.genfromtxt(frame_as_list_text_for_bonds)
            lammps_dix['index'][frames_cnt]['bonds'] = bond_locations
            # print(bond_locations)
            frames_cnt += 1

    return lammps_dix


if __name__ == "__main__":

    # fname = "propagation.lammpstrj"
    fname = "propagation_small.lammpstrj"

    lammps_dix = load_lammps(fname)
