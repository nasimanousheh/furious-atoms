from furiousatoms.parsers.pdb_parser import PDBParser #TODO remove
from furiousatoms.molecular import MolecularStructure
from furiousatoms.parsers.parser_util import float_or_zero
from furiousatoms.savers.base_saver import BaseSaver
import numpy as np

ATOM_FORMAT = "{label:<6s}{serial:5d}  {name:<4s}{altLoc:1s}{resName:<3s}" + \
    "{chainID:1s}{resSeq:>4s}{iCode:1s}" + \
    "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}" + \
    "{tempFactor:6.2f}      {segID:<4s}{element:>2s}{charge:2s}\n"

DEFAULT_HEADER = "REMARK Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"
DEFAULT_LABEL = "ATOM"
DEFAULT_RES_NAME = "MOL"
DEFAULT_ALT_LOC = " "
DEFAULT_CHAIN_ID = "X"
DEFAULT_RES_SEQ = "1"
DEFAULT_ICODE = "1"
DEFAULT_OCCUPANCY = 1.0
DEFAULT_TEMP_FACTOR = 0.0
DEFAULT_SEG_ID = "    "
DEFAULT_CHARGE = "  "

class PDBSaver(BaseSaver):
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
    
    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        http://www.bmsc.washington.edu/CrystaLinks/man/pdb/

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.saved_atom_ids = np.zeros(shape=(len(structure.pos)))
        self.atom_serial = 1
        header = old_fp.readline()
        new_fp.write(header)
        
        while True:
            old_line = old_fp.readline()
            if not old_line:
                break
            line_upper = old_line.upper()
            length = len(old_line)
            if line_upper.startswith("CONECT"):
                pass
            elif line_upper.startswith("ATOM") or line_upper.startswith("HETATM"):
                if length >= 54: #does it contain position?
                    atom_id = int(old_line[7:11]) - 1
                    self.saved_atom_ids[atom_id - 1] = 1

                    #Save each field so that we can use them to "guess" 
                    # fields for newly-added atoms later.
                    self.last_label = old_line[:6]
                    self.last_altLoc = old_line[17]
                    self.last_resName = old_line[18:20]
                    self.last_chainID = old_line[21:22]
                    self.last_resSeq = old_line[23:26]
                    self.last_iCode = old_line[27]
                    #Below fields may not be present
                    self.last_occupancy = float_or_zero(old_line[55:60]) if length >= 60 else DEFAULT_OCCUPANCY
                    self.last_tempFactor = float_or_zero(old_line[61:66]) if length >= 66 else DEFAULT_TEMP_FACTOR
                    self.last_segID = old_line[73:76] if length >= 76 else DEFAULT_SEG_ID
                    self.last_charge = old_line[79:80] if length >= 80 else DEFAULT_CHARGE
                    
                    self._write_atom(new_fp, structure, atom_id)
            elif line_upper.startswith("CRYST1"):
                #Write box size in angstroms but keep rest of old line
                box_size_format = "CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f} "
                new_fp.write(box_size_format.format(
                    box=structure.box_size,
                ))
                if len(old_line) > 34:
                    new_fp.write(old_line[34:])
                else:
                    new_fp.write('\n')
            elif not line_upper.startswith("END"): #keep REMARKs and other lines
                new_fp.write(old_line)

        #Write new atoms
        self.last_label = DEFAULT_LABEL
        for atom_id in range(len(structure.pos)):
            if not self.saved_atom_ids[atom_id]:
                self._write_atom(new_fp, structure, atom_id)

        self._write_all_bonds(new_fp, structure)

        new_fp.write("END\n")


    def save_to_file_use_defaults(self, fname, structure: MolecularStructure):
        self.atom_serial = 1
        self.last_label=DEFAULT_LABEL
        self.last_altLoc=DEFAULT_ALT_LOC
        self.last_resName=DEFAULT_RES_NAME
        self.last_chainID=DEFAULT_CHAIN_ID
        self.last_resSeq=DEFAULT_RES_SEQ
        self.last_iCode=DEFAULT_ICODE
        self.last_occupancy=DEFAULT_OCCUPANCY
        self.last_tempFactor=DEFAULT_TEMP_FACTOR
        self.last_segID=DEFAULT_SEG_ID
        self.last_charge=DEFAULT_CHARGE

        with open(fname, 'w') as fp:
            fp.write(DEFAULT_HEADER)

            #Write box size
            boxSizeFormat = "CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}" + \
                   "{ang[0]:7.2f}{ang[1]:7.2f}{ang[2]:7.2f} " + \
                   "{spacegroup:<11s}{zvalue:4d}\n"
            fp.write(boxSizeFormat.format(
                box=structure.box_size, #size in angstroms
                ang=[0, 0, 0], #alpha/beta/gamma (degrees)
                spacegroup="P",
                zvalue=1
            ))

            for i in range(len(structure.pos)):
                self._write_atom(fp, structure, i)

            self._write_all_bonds(fp, structure)

            fp.write("END\n")


    def _write_atom(self, new_fp, structure, atom_id):
        if self.is_atom_deleted(atom_id):
            return
        new_fp.write(ATOM_FORMAT.format(
            label=self.last_label,
            serial=self.atom_serial,
            name=structure.atom_types[atom_id],
            altLoc=self.last_altLoc,
            resName=self.last_resName,
            chainID=self.last_chainID,
            resSeq=self.last_resSeq,
            iCode=self.last_iCode,
            pos=structure.pos[atom_id],
            occupancy=self.last_occupancy,
            tempFactor=self.last_tempFactor,
            segID=self.last_segID,
            element=structure.atom_types[atom_id],
            charge=self.last_charge
        ))
        self.atom_serial += 1


    def _write_all_bonds(self, new_fp, structure):
        bond_map = {}
        i = 0
        for bond in self.corrected_bonds:
            if not self.is_bond_deleted(i):
                bond_list = bond_map.get(bond[0] + 1, [])
                bond_list.append(bond[1] + 1)
                bond_map[bond[0] + 1] = bond_list
                #Assume all bonds go in both directions
                bond_list = bond_map.get(bond[1] + 1, [])
                bond_list.append(bond[0] + 1)
                bond_map[bond[1] + 1] = bond_list
            i += 1

        for atom_id in range(1, len(structure.pos) + 1):
            bond_list = bond_map.get(atom_id, [])
            if len(bond_list) > 0:
                conect = ["{0:5d}".format(entry) for entry in bond_list]
                conect = "".join(conect)
                bond_format = "CONECT{atomId:5d}{conect}\n"
                new_fp.write(bond_format.format(atomId=atom_id, conect=conect))


#TODO remove
if __name__ == "__main__":
    # INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\furiousatoms\\tests\\test_data\\CB_18\\CB_18.pdb"
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\test.pdb"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\test.pdb"
        
    structure = PDBParser().parse(INPUT_FNAME)
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_particles[5] = True
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = PDBSaver(deleted_particles, deleted_bonds)
    saver.save_to_file(OUTPUT_FNAME, INPUT_FNAME, structure)
    # saver.save_to_file_use_defaults(OUTPUT_FNAME, structure)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')