# Standard package
import os
import sys

import furiousatoms
from PySide2 import QtCore, QtUiTools
from warnings import warn

from fury.lib import Texture, ImageReader2Factory, ImageFlip
from furiousatoms.molecular import MolecularStructure
from furiousatoms.parsers.gromacs_parser import GROMACSParser
from furiousatoms.parsers.lammps_parser import LAMMPSParser
from furiousatoms.parsers.pdb_parser import PDBParser
from furiousatoms.parsers.xyz_parser import XYZParser
from furiousatoms.savers.gromacs_saver import GROMACSSaver
from furiousatoms.savers.lammps_saver import LAMMPSSaver
from furiousatoms.savers.pdb_saver import PDBSaver
from furiousatoms.savers.xyz_saver import XYZSaver


def is_frozen():
    return getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')


def get_frozen_path():
    base_path = ''
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return base_path


def get_package_path():
    """Returns package path."""
    package_path = os.path.dirname(os.path.realpath(furiousatoms.__file__))
    # print("package path : {0}".format(package_path))
    return package_path


def get_application_path():
    """Returns application path."""
    app_path = get_package_path()
    if is_frozen():
        app_path = os.path.dirname(sys.executable)
    return app_path


def get_languages_path():
    base_path = get_frozen_path() if is_frozen() else get_application_path()
    l_path = os.path.join(base_path, "languages")
    if not os.path.isdir(l_path):
        os.mkdir(l_path)
    return l_path


def get_resources_file(fname):
    base_path = get_frozen_path() if is_frozen() else get_application_path()
    return os.path.join(base_path, "resources", fname)


def load_ui_widget(uifilename, cls_to_register=None, parent=None):
    loader = QtUiTools.QUiLoader()
    loader.setLanguageChangeEnabled(True)
    loader.setTranslationEnabled(True)
    if cls_to_register:
        for cls in cls_to_register:
            loader.registerCustomWidget(cls)

    base_path = get_frozen_path() if is_frozen() else get_application_path()
    ui_fpath = os.path.join(base_path, "forms", uifilename)
    if not os.path.isfile(ui_fpath):
        msg = '{} does not exist'.format(ui_fpath)
        raise ValueError(msg)

    print("Loading UI file: {0}".format(os.path.join(base_path,
                                        "forms", uifilename)))
    uifile = QtCore.QFile(ui_fpath)
    uifile.open(QtCore.QFile.ReadOnly | QtCore.QIODevice.Text)
    uifile.reset()
    ui = loader.load(uifile, parent)
    uifile.close()
    return ui


def load_files(fname):
    #Choose parser based on file extension
    if fname.endswith(".pdb") or fname.endswith(".cc1"):
        return PDBParser().parse(fname)
    elif fname.endswith(".data") or fname.endswith(".dat") or fname.endswith(".lmp"):
        return LAMMPSParser().parse(fname)
    elif fname.endswith(".xyz"):
        return XYZParser().parse(fname)
    elif fname.endswith(".gro"):
        return GROMACSParser().parse(fname)
    
    #Default: try lots of parsers and guess the format
    best_structure  = MolecularStructure.create_empty()
    for parser in (PDBParser(), LAMMPSParser(), XYZParser(), GROMACSParser()):
        structure = parser.parse(fname)
        if len(structure.bonds) > len(best_structure .bonds) or 0 == len(best_structure .pos):
            best_structure  = structure
    return best_structure 

def save_file(new_path, old_path, structure, deleted_particles, deleted_bonds):
    saver = None
    if old_path.endswith(".pdb") or old_path.endswith(".cc1"):
        saver = PDBSaver(deleted_particles, deleted_bonds)
    elif old_path.endswith(".data") or old_path.endswith(".dat") or old_path.endswith(".lmp"):
        saver  = LAMMPSSaver(deleted_particles, deleted_bonds)
    elif old_path.endswith(".xyz"):
        saver = XYZSaver(deleted_particles, deleted_bonds)
    elif old_path.endswith(".gro"):
        saver = GROMACSSaver(deleted_particles, deleted_bonds)
    
    if not saver:
        warn("Could not decide on which saver to use. Using the default format (XYZ).")
        saver = XYZSaver(deleted_particles, deleted_bonds)
        
    saver.save_to_file(new_path, old_path, structure)


def read_cubemap(folderRoot, fileRoot, ext, key):
    """Read the cube map.

    Parameters
    ----------
    folderRoot: str
        The folder where the cube maps are stored.
    fileRoot : str
        The root of the individual cube map file names.
    ext : str
        The extension of the cube map files.
    key: str
        The key to data used to build the full file name.

    Returns
    -------
    texture : Texture
        The cubemap texture.

    """
    # A map of cube map naming conventions and the corresponding file name
    # components.
    fileNames = {
        0: ['right', 'left', 'top', 'bottom', 'front', 'back'],
        1: ['posx', 'negx', 'posy', 'negy', 'posz', 'negz'],
        2: ['-px', '-nx', '-py', '-ny', '-pz', '-nz'],
        3: ['0', '1', '2', '3', '4', '5']}
    if key in fileNames:
        fns = fileNames[key]
    else:
        print('ReadCubeMap(): invalid key, unable to continue.')
        sys.exit()
    texture = Texture()
    texture.CubeMapOn()
    # Build the file names.
    for i in range(0, len(fns)):
        fns[i] = folderRoot + fileRoot + fns[i] + ext
        if not os.path.isfile(fns[i]):
            print('Nonexistent texture file:', fns[i])
            return texture
    i = 0
    for fn in fns:
        # Read the images
        readerFactory = ImageReader2Factory()
        imgReader = readerFactory.CreateImageReader2(fn)
        imgReader.SetFileName(fn)

        flip = ImageFlip
        flip.SetInputConnection(imgReader.GetOutputPort())
        flip.SetFilteredAxis(1)  # flip y axis
        texture.SetInputConnection(i, flip.GetOutputPort(0))
        i += 1
    return texture
