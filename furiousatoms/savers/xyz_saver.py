import numpy as np

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class XYZSaver():
    def write_all_lines(self, fname, box_size: list, positions: np.array, bonds: np.array, atom_types: np.array):
        with open(fname, 'w') as fp:
            fp.write("%d\n"%positions.shape[0])
            fp.write(HEADER)

            for i in range(len(positions)):
                x, y, z = positions[i]
                fp.write("%2s %8.6f %8.6f %8.6f\n"%(atom_types[i], x, y, z))