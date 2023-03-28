import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero, has_bond
from furiousatoms.element_lookup import lookup_symbol_by_mass
from warnings import warn


class LAMMPSParser(BaseParser):
    def __init__(self) -> None:
        super().__init__()
        self.id_to_symbol_dict = {}
    
    def parse_mass(self, line):
        words = line.split()
        key = words[0]
        mass = float(words[1])
        try:
            self.id_to_symbol_dict[key] = lookup_symbol_by_mass(mass)
        except ValueError as ex:
            self.id_to_symbol_dict[key] = mass
    
    def get_atom_type(self, typeKey):
        '''
        This method accounts for the fact that some atoms are identified by integers 
        rather than their symbol. For example, given the masses:
        ```
        1 1.0   (Hydrogen)
        2 12.0   (Carbon)
        ```
        The file's atoms might be:
        ```
        1 2 -5.46 0.0 0.0
        2 1 -3.9 0.0 0.0
        3 2 -3.12 -1.351 0.0
        ```
        Parameter `typeKey` can either be the `int` that identifies the element
        or the element symbol itself. In either case, `typeKey` must have been
        read in by parse_mass(...) first.
        '''
        if typeKey in self.id_to_symbol_dict:
            return self.id_to_symbol_dict[typeKey]
        return typeKey

    def parse_atom(self, line):
        '''
        Atom lines come in many styles, but we assume the XYZ coordinates are the last 3 
        values as per the 'atomic' and 'full' styles. 
        '''
        try:
            pos = np.zeros((3))
            typ = ''
            words = line.split()
            if len(words) >= 7: #'Full' style
                typ = self.get_atom_type(words[2])
                pos[0] = float_or_zero(words[4])
                pos[1] = float_or_zero(words[5])
                pos[2] = float_or_zero(words[6])
            elif len(words) >= 5: #'Atomic' style
                typ = self.get_atom_type(words[1])
                pos[0] = float_or_zero(words[2])
                pos[1] = float_or_zero(words[3])
                pos[2] = float_or_zero(words[4])
            else:
                raise IndexError("Line too short")
            self.atom_types.append(typ)
            self.positions.append(pos)
        except IndexError:
            self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.line_id)
        except:
            self.errors += "Failure processing line #%d.\n"%(self.line_id)
    
    def parse_bond(self, line):
        words = line.split()
        try:
            if len(words) >= 4:
                bond = np.zeros((2), dtype='int')
                #Record IDs of the two atoms to connect.
                bond[0] = int(words[2]) - 1
                bond[1] = int(words[3]) - 1
                if bond[0] < 0 or bond[1] < 0:
                    raise ValueError("Atom ID cannot be zero or negative.")
                if bond[0] != bond[1] and not has_bond(self.bonds, bond):
                    self.bonds.append(bond)
            else: 
                raise IndexError("Line too short")
        except:
            if len(words) < 2:
                self.errors += "Line #%d is too short. At least two atom IDs are needed for a bond.\n"%(self.line_id)
            else:
                self.errors += "Could not read a bond from line #%d since it is not well-formatted. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.line_id)

    def parse_box_size(self, line):
        '''
        A box size dimension will look like this:
        ` -2.365000 7.095000  xlo xhi` 
        This function checks what label is at the end of the line and extracts the hi dimension.
        '''
        try:
            words = line.split()
            trimmed = line.strip()
            if len(words) >= 4: 
                #Always parse the "hi" dimension
                if trimmed.endswith("xlo xhi"):
                    self.box_size[0] = float_or_zero(words[1])
                elif trimmed.endswith("ylo yhi"):
                    self.box_size[1] = float_or_zero(words[1])
                elif trimmed.endswith("zlo zhi"):
                    self.box_size[2] = float_or_zero(words[1])
                else:
                    raise ValueError("Invalid format")
            else:
                self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.line_id)
        except:
            self.errors += "Failure processing line #%d.\n"%(self.line_id)

    def skip_line(self, line):
        pass

    def parse_line(self, line):
        #Choose how to parse subsequent lines
        line = line.strip()
        header = line.lower()
        if len(header) == 0:
            pass
        elif header.startswith("masses"):
            self.parser_method = self.parse_mass
        elif header.startswith("atoms"):
            self.parser_method = self.parse_atom
        elif header.startswith("bonds"):
            self.parser_method = self.parse_bond
        elif header.endswith("hi") and "lo" in header:
            #Special case: box size doesn't use a header
            self.parse_box_size(line)
        elif "coeffs" in header:
            self.parser_method = self.skip_line
        elif len(header.split()) == 1:
            #For non-implemented headers like Angles and Dihedrals
            self.parser_method = self.skip_line
        
        #Perform the parse
        elif len(header.strip()) > 0:
            if self.parser_method != None:
                self.parser_method(line)
            else:
                self.errors += "Unable to process line #%d.\n"%(self.line_id)
        self.line_id += 1