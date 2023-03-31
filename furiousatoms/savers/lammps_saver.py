import numpy as np
from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver
from furiousatoms.parsers.parser_util import float_or_zero
from warnings import warn

ATOM_FORMAT = "%d %2s %f %f %f\n"
BOND_FORMAT = "%d %s %d %d\n"
DEFAULT_HEADER = "LAMMPS data file. Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"
DEFAULT_BOND_TYPE = 1

class LAMMPSSaver(BaseSaver):
    def write_single_atom(self, fp, structure, atom_id):
        if not self.is_atom_deleted(atom_id):
            x, y, z = structure.pos[atom_id]
            typ = structure.atom_types[atom_id]
            fp.write(ATOM_FORMAT%(self.atom_serial, typ, x, y, z))
            self.atom_serial += 1
            self.saved_atom_ids[atom_id] = 1

    def write_single_bond(self, fp, bond_id):
        if not self.is_bond_deleted(bond_id):
            atom1, atom2 = self.corrected_bonds[bond_id]
            atom1 += 1
            atom2 += 1
            fp.write(BOND_FORMAT%(self.bond_serial, self.bond_type, atom1, atom2))
            self.bond_serial += 1

    def write_all_atoms(self, new_fp, old_fp, structure):
        '''
        Always write in 'atomic' style.
        Syntax: `atom-ID atom-type x y z`
        '''
        self.saved_atom_ids = np.zeros(shape=(len(structure.pos)))
        i = 0
        while True:
            prev_pos = old_fp.tell()
            line = old_fp.readline()
            words = line.split()
            if len(words) == 0:
                continue
            if len(words) < 5: #num. fields in 'atomic' style
                old_fp.seek(prev_pos)
                break
            atom_id = int(words[0]) - 1
            self.write_single_atom(new_fp, structure, atom_id)
            i += 1

        #Write any atoms missing from the source file
        for atom_id in range(len(structure.pos)):
            if not self.saved_atom_ids[atom_id]:
                self.write_single_atom(new_fp, structure, atom_id)

        new_fp.write("\n")


    def write_all_bonds(self, new_fp, old_fp):
        i = 0
        self.bond_serial = 1
        while True:
            prev_pos = old_fp.tell()
            line = old_fp.readline()
            words = line.split()
            if len(words) == 0:
                continue
            if len(words) < 4 or i > len(self.corrected_bonds):
                old_fp.seek(prev_pos)
                break
            
            self.bond_type = words[1]
            self.write_single_bond(new_fp, i)
            i += 1
        new_fp.write("\n")
        

    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://www.smcm.iqfr.csic.es/docs/lammps/read_data.html

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.saved_atom_ids = np.zeros(shape=(len(structure.pos)))
        self.atom_serial = 1
        should_keep_line = True
        new_fp.write(old_fp.readline()) #file comment
        
        while True:
            old_line = old_fp.readline()
            if not old_line:
                break
            line = old_line.strip()
            line_lower = line.lower()

            #LAMMPS BODY ##############################
            if line_lower.startswith("masses"): 
                should_keep_line = True
            elif line_lower.startswith("atoms"): 
                new_fp.write("Atoms          #atomic\n\n")
                self.write_all_atoms(new_fp, old_fp, structure)
                should_keep_line = False
            elif line_lower.startswith("bonds"): 
                new_fp.write("Bonds\n\n")
                self.write_all_bonds(new_fp, old_fp)
                should_keep_line = False
            #LAMMPS HEADER ##############################
            elif line_lower.endswith("atoms"):
                num_atoms = len(np.argwhere(self.deleted_particles == False))
                new_fp.write("%d\tatoms\n"%(num_atoms))
                should_keep_line = False
            elif line_lower.endswith("bonds"):
                num_bonds = len(np.argwhere(self.deleted_bonds == False))
                new_fp.write("%d\tbonds\n"%(num_bonds))
                should_keep_line = False
            elif line_lower.endswith("atom types"):
                num_types = len(np.unique(structure.atom_types))
                new_fp.write("%d\tatom types\n"%(num_types))
                should_keep_line = False
            elif line_lower.endswith("hi") and "lo" in line_lower:
                words = line.split()
                if len(words) != 4:
                    new_fp.write(line) #handles invalid box size lines
                    continue
                lo = float_or_zero(words[0])
                if line_lower.endswith("xlo xhi"):
                    hi = structure.box_size[0]
                elif line_lower.endswith("ylo yhi"):
                    hi = structure.box_size[1]
                elif line_lower.endswith("zlo zhi"):
                    hi = structure.box_size[2]
                else:
                    hi = 0.0
                new_fp.write("%f\t%f\t%s %s\n"%(lo, hi, words[2], words[3]))
                should_keep_line = False
            else:
                should_keep_line = True

            
            if should_keep_line:
                new_fp.write(line)
                if "\n" not in line:
                    new_fp.write("\n")


    def _save_to_file_use_defaults(self, new_fp, structure: MolecularStructure):
        from furiousatoms.element_lookup import lookup_mass_by_symbol, DEFAULT_MASS
        '''
        Write all data to a file using the format specified in 
        https://www.smcm.iqfr.csic.es/docs/lammps/read_data.html

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.saved_atom_ids = np.zeros(shape=(len(structure.pos)))
        self.atom_serial = 1
        new_fp.write(DEFAULT_HEADER + "\n")

        new_fp.write("%d\tatoms\n"%(len(structure.pos)))
        new_fp.write("%d\tbonds\n"%(len(structure.bonds)))
        num_types = len(np.unique(structure.atom_types))
        new_fp.write("%d\tatom types\n"%(num_types))
        DEFAULT_NUM_BOND_TYPES = 1
        new_fp.write("%d\tbond types\n"%(DEFAULT_NUM_BOND_TYPES))
        new_fp.write("\n")

        for i, label in enumerate(("xlo xhi", "ylo yhi", "zlo zhi")):
            lo = structure.box_size[i]
            hi = lo
            new_fp.write("%f\t%f\t%s\n"%(lo, hi, label))

        new_fp.write("\nMasses\n\n")
        for typ in np.unique(structure.atom_types):
            try:
                mass = lookup_mass_by_symbol(typ)
                new_fp.write("%s %f\n"%(typ, mass))
            except ValueError as ex:
                warn(str(ex))
                new_fp.write("%s %f\n"%(typ, DEFAULT_MASS))
        
        new_fp.write("\nAtoms          #atomic\n\n")
        for atom_id in range(len(structure.pos)):
            self.write_single_atom(new_fp, structure, atom_id)

        new_fp.write("\nBonds\n\n")
        self.bond_serial = 1
        self.bond_type = DEFAULT_BOND_TYPE
        for bond_id in range(len(self.corrected_bonds)):
            self.write_single_bond(new_fp, bond_id)