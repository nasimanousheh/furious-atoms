import numpy as np
from furiousatoms.parsers.gromacs_parser import GROMACSParser #TODO remove

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class GROMACSSaver():
    '''
    Write all data to a file using the format specified in 
    https://manual.gromacs.org/archive/5.0.3/online/gro.html
    '''
    def writeAllLines(self, fname, boxSize: list, positions: np.array, bonds: np.array, atomTypes: np.array):
        with open(fname, 'w') as fp:
            fp.write(HEADER)
            fp.write("%d\n"%positions.shape[0]) #atom count

            for i in range(len(positions)):
                x, y, z = positions[i]
                defaultResidueNumber = 2
                defaultResidueName = "MOL"
                atomNumber = 0 #TODO add VTK table lookup
                fp.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(
                    defaultResidueNumber,
                    defaultResidueName,
                    atomTypes[i], 
                    atomNumber,
                    x, y, z
                ))
            
            fp.write("%4.5f %4.5f %4.5f\n"%(boxSize[0], boxSize[1], boxSize[2]))

#TODO remove
if __name__ == "__main__":
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\Graphdiyne\\graphdiyne_unitcell.gro"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\graphdiyne_test.gro"
    boxSize, positions, bonds, atomTypes = GROMACSParser().parse(INPUT_FNAME)
    saver = GROMACSSaver()
    saver.writeAllLines(OUTPUT_FNAME, boxSize, positions, bonds, atomTypes)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')