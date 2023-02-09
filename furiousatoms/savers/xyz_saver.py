import numpy as np
from furiousatoms.parsers.xyz_parser import XYZParser #TODO remove

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class XYZSaver():
    def writeAllLines(self, fname, boxSize: list, positions: np.array, bonds: np.array, atomTypes: np.array):
        with open(fname, 'w') as fp:
            fp.write("%d\n"%positions.shape[0])
            fp.write(HEADER)

            for i in range(len(positions)):
                x, y, z = positions[i]
                fp.write("%2s %8.6f %8.6f %8.6f\n"%(atomTypes[i], x, y, z))


#TODO remove
if __name__ == "__main__":
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\Graphdiyne\\graphdiyne_unitcell.xyz"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\graphdiyne_test.xyz"
    boxSize, positions, bonds, atomTypes = XYZParser().parse(INPUT_FNAME)
    saver = XYZSaver()
    saver.writeAllLines(OUTPUT_FNAME, boxSize, positions, bonds, atomTypes)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')