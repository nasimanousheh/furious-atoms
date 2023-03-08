import numpy as np

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class GROMACSSaver():
    '''
    Write all data to a file using the format specified in 
    https://manual.gromacs.org/archive/5.0.3/online/gro.html
    '''
    def write_all_lines(self, fname, box_size: list, positions: np.array, bonds: np.array, atom_types: np.array):
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
                    atom_types[i], 
                    atomNumber,
                    x, y, z
                ))
            
            fp.write("%4.5f %4.5f %4.5f\n"%(box_size[0], box_size[1], box_size[2]))