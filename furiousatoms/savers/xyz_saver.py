import numpy as np

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class XYZSaver():
    def writeAllLines(self, fname, boxSize: list, positions: np.array, bonds: np.array, atomTypes: np.array):
        with open(fname, 'w') as fp:
            fp.write("%d\n"%positions.shape[0])
            fp.write(HEADER)

            for i in range(len(positions)):
                x, y, z = positions[i]
                fp.write("%2s %8.6f %8.6f %8.6f\n"%(atomTypes[i], x, y, z))