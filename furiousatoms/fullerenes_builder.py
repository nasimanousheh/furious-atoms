import numpy as np


def load_CC1_file(fname, debug=False):
    lammps_file = open(fname, 'r')
    lines = lammps_file.readlines()
    no_lines = len(lines)
    lammps_dix = {'index': {}}
    frames_cnt = 0
    no_atoms = int(len(lines)-1)
    for i in range(no_lines-1):
        lammps_dix['index'][frames_cnt] = {'no_atoms': no_atoms, 'coords': [], 'bonds': []}
        frame_as_list_text = lines[i+1: i+1 + no_atoms]
        frame_positions = np.genfromtxt(frame_as_list_text)
        lammps_dix['index'][frames_cnt]['coords'] = frame_positions
        frames_cnt += 1

    frames_cnt = 0
    no_atoms = lammps_dix['index'][frames_cnt]['no_atoms']
    pos = lammps_dix['index'][frames_cnt]['coords'][:, 2:5]
    bonds = np.asarray(lammps_dix['index'][frames_cnt]['coords'][:, 6:9], dtype="i8")
    file_name = 'fullerene.pdb'
    outdump = open(file_name, "w")
    outdump.write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P1\n")
    for i in range(no_atoms):
        tmp = "{:6s}{:5d}  {:5s}{:2s}{:3s} {:2s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM', 1 + i, 'C', 'C', 'A', '1', float("{:.3f}".format(pos[i][0])), float("{:.3f}".format(pos[i][1])),float("{:.3f}".format(pos[i][2])),1,0)
        outdump.write(tmp)
    outdump.write("TER\n")
    for j in range(no_atoms):
        tmp_bond = "{:4s}{:5d} {:4d} {:4d} {:4d}\n".format('CONECT', 1 + j, bonds[j][0], bonds[j][1], bonds[j][2])
        outdump.write(tmp_bond)
    outdump.write("{:4s} {:5d} {:5d} {:5d} {:5d}\n".format('MASTER        0    0    0    0    0    0    0    0', no_atoms, 0, no_atoms, 0))
    outdump.write("END\n")
    return file_name
