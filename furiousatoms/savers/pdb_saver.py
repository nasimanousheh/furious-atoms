import numpy as np
from furiousatoms.parsers.pdb_parser import PDBParser #TODO remove

HEADER = "REMARK Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class PDBSaver():
    '''
    .. Format strings and the following links were taken MDAnalysis.
    .. https://github.com/MDAnalysis/mdanalysis/ 
    
    .. _`PDB 3.3 standard`:
       http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    .. _ATOM: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    .. _COMPND: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#COMPND
    .. _CONECT: http://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT
    .. _END: http://www.wwpdb.org/documentation/file-format-content/format33/sect11.html#END
    .. _ENDMDL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ENDMDL
    .. _HEADER: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#HEADER
    .. _HETATM: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
    .. _MODEL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL
    .. _NUMMDL: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#NUMMDL
    .. _REMARKS: http://www.wwpdb.org/documentation/file-format-content/format33/remarks.html
    .. _TITLE: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#TITLE
    '''
    
    def writeAllLines(self, fname, boxSize: list, positions: np.array, bonds: np.array, atomTypes: np.array):
        with open(fname, 'w') as fp:
            fp.write(HEADER)

            #Write box size
            #FIXME: CRYST1 also needs alpha/beta/gamma (degrees), space group (string), and z (integer)
            boxSizeFormat = "CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}" + \
                   "{ang[0]:7.2f}{ang[1]:7.2f}{ang[2]:7.2f} " + \
                   "{spacegroup:<11s}{zvalue:4d}\n"
            fp.write(boxSizeFormat.format(
                box=boxSize, #size in angstroms
                ang=[0, 0, 0], #alpha/beta/gamma (degrees)
                spacegroup="P",
                zvalue=1
            ))

            #Write atoms
            for i in range(len(positions)):
                atomFormat = "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}" + \
                    "{chainID:1s}{resSeq:4d}{iCode:1s}" + \
                    "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}" + \
                    "{tempFactor:6.2f}      {segID:<4s}{element:>2s}{charge:2s}\n"
                fp.write(atomFormat.format(
                    serial=i,
                    name=atomTypes[i],
                    altLoc=" ",
                    resName=" MOL",
                    chainID=" ",
                    resSeq=1,
                    iCode="1",
                    pos=positions[i],
                    occupancy=1.00,
                    tempFactor=0,
                    segID="    ",
                    element=atomTypes[i],
                    charge="  "
                ))

            #Organize bonds
            bondMap = {}
            for bond in bonds:
                bondList = bondMap.get(bond[0] + 1, [])
                bondList.append(bond[1] + 1)
                bondMap[bond[0] + 1] = bondList

            #Write bonds
            for atomId, bondList in bondMap.items():
                conect = ["{0:5d}".format(entry) for entry in bondList]
                conect = "".join(conect)
                bondFormat = "CONECT{atomId:5d}{conect}\n"
                fp.write(bondFormat.format(atomId=atomId, conect=conect))

            fp.write("END\n")


#TODO remove
if __name__ == "__main__":
    # INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\Graphdiyne\\graphdiyne_unitcell.pdb"
    # OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\graphdiyne_test.pdb"
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\CB_18.pdb"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\examples\\CB_18_test.pdb"
    boxSize, positions, bonds, atomTypes = PDBParser().parse(INPUT_FNAME)
    saver = PDBSaver()
    saver.writeAllLines(OUTPUT_FNAME, boxSize, positions, bonds, atomTypes)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')