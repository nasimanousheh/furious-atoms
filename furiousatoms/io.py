# Standard package
import os
import sys

import furiousatoms
from PySide2 import QtCore, QtUiTools
import MDAnalysis

from fury.lib import Texture, ImageReader2Factory, ImageFlip


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


def get_dateset_file(relative_path):
    base_path = ''
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    # base_path = get_frozen_path() if is_frozen() else get_application_path()
    return os.path.join(base_path, relative_path)

# def get_dateset_file():
#     base_path = get_frozen_path() if is_frozen() else get_application_path()
#     # return os.path.join(base_path, "fullerene_dataset")
#     dataset_path = os.path.join(base_path, "fullerene_dataset")
#     if not os.path.isdir(dataset_path):
#         os.mkdir(dataset_path)
#     return dataset_path

# def get_dateset_file():
#     base_path = get_frozen_path() if is_frozen() else get_application_path()
#     l_path = os.path.join(base_path, "fullerene_dataset")
#     if not os.path.isdir(l_path):
#         os.mkdir(l_path)
#     return l_path


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


def load_files(fname, debug=False):
    load_file_export = open(fname, 'r')
    lines = load_file_export.readlines()
    no_lines = len(lines)
    frames_cnt = 0
    format_data = None
    no_bonds = 0
    bonds = 0
    for i in range(no_lines):
        line = lines[i]
        if 'item: number of atoms'.upper() in line:
            format_data = 'LAMMPSDUMP'
            break
        i += 1

    load_file = MDAnalysis.Universe(fname, format=format_data)
    if format_data is None:
        load_file = MDAnalysis.Universe(fname)
        try:
            bonds = load_file.bonds.to_indices()
            no_bonds = len(load_file.bonds)
        except MDAnalysis.exceptions.NoDataError:
            no_bonds = 0
    return load_file, no_bonds


def create_universe(pos, bonds, atom_types, box_lx, box_ly, box_lz):
    num_atoms = pos.shape[0]
    universe = MDAnalysis.Universe.empty(num_atoms, trajectory=True, n_residues=1)
    universe.atoms.positions = pos
    n_residues = 1
    atom_types_list = list(atom_types)
    universe.add_TopologyAttr('name', atom_types_list)
    universe.add_TopologyAttr('type', atom_types_list)
    universe.add_TopologyAttr('resname', ['MOL']*n_residues)
    universe.add_TopologyAttr('masses')
    universe.trajectory.ts.dimensions = [box_lx, box_ly, box_lz, 90, 90, 90]
    try:
        universe.add_bonds(bonds)
    except:
        pass
    cog = universe.atoms.center_of_geometry()
    universe.atoms.positions -= cog
    return universe

def merged_two_universes(pos_uni_1, bonds_uni_1, atom_types_uni_1, pos_uni_2, bonds_uni_2, atom_types_uni_2, box_lx, box_ly, box_lz):
    num_atoms_1 = pos_uni_1.shape[0]
    universe_1 = MDAnalysis.Universe.empty(num_atoms_1, trajectory=True, n_residues=1)
    universe_1.atoms.positions = pos_uni_1
    n_residues = 1
    atom_types_list_1 = list(atom_types_uni_1)
    universe_1.add_TopologyAttr('name', atom_types_list_1)
    universe_1.add_TopologyAttr('type', atom_types_list_1)
    universe_1.add_TopologyAttr('resname', ['MOL']*n_residues)
    universe_1.add_bonds(bonds_uni_1)
    num_atoms_2 = pos_uni_2.shape[0]
    universe_2 = MDAnalysis.Universe.empty(num_atoms_2, trajectory=True, n_residues=1)
    universe_2.atoms.positions = pos_uni_2
    universe_2.add_bonds(bonds_uni_2)
    n_residues = 1
    universe_2.add_TopologyAttr('name', atom_types_uni_2)
    universe_2.add_TopologyAttr('type', atom_types_uni_2)
    universe_2.add_TopologyAttr('resname', ['MOL']*n_residues)
    merged_universes = MDAnalysis.Merge(universe_1.atoms, universe_2.atoms)
    merged_universes.trajectory.ts.dimensions = [box_lx, box_ly, box_lz, 90, 90, 90]
    merged_universes.add_bonds(universe_1.bonds.indices)
    return merged_universes

def merged_universe_with_H(pos_uni_1, bonds_uni_1, atom_types_uni_1, pos_uni_2, bonds_uni_2, atom_types_uni_2, box_lx, box_ly, box_lz):
    num_atoms_1 = pos_uni_1.shape[0]
    universe_1 = MDAnalysis.Universe.empty(num_atoms_1, trajectory=True, n_residues=1)
    universe_1.atoms.positions = pos_uni_1
    n_residues = 1
    atom_types_list_1 = list(atom_types_uni_1)
    universe_1.add_TopologyAttr('name', atom_types_list_1)
    universe_1.add_TopologyAttr('type', atom_types_list_1)
    universe_1.add_TopologyAttr('resname', ['MOL']*n_residues)
    universe_1.add_bonds(bonds_uni_1)
    num_atoms_2 = pos_uni_2.shape[0]
    universe_2 = MDAnalysis.Universe.empty(num_atoms_2, trajectory=True, n_residues=1)
    universe_2.atoms.positions = pos_uni_2
    # universe_2.add_bonds(bonds_uni_2)
    n_residues = 1
    universe_2.add_TopologyAttr('name', atom_types_uni_2)
    universe_2.add_TopologyAttr('type', atom_types_uni_2)
    universe_2.add_TopologyAttr('resname', ['MOL']*n_residues)
    merged_universe_Hyd = MDAnalysis.Merge(universe_1.atoms, universe_2.atoms)
    merged_universe_Hyd.add_bonds(bonds_uni_1)
    merged_universe_Hyd.add_bonds(bonds_uni_2)
    try:
        box_lx or box_ly or box_lz
    except NameError:
        box_lx = box_ly = box_lz = 0.0
    merged_universe_Hyd.trajectory.ts.dimensions = [box_lx, box_ly, box_lx, 90, 90, 90]
    cog = merged_universe_Hyd.atoms.center_of_geometry()
    merged_universe_Hyd.atoms.positions -= cog
    return merged_universe_Hyd
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
