# Standard package
import os
import sys

import furiousatoms

from PySide2 import QtCore, QtUiTools
import MDAnalysis


def get_package_path():
    """Returns package path."""
    package_path = os.path.dirname(os.path.realpath(furiousatoms.__file__))
    # print("package path : {0}".format(package_path))
    return package_path


def get_application_path():
    """Returns application path."""
    app_path = get_package_path()
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        app_path = os.path.dirname(sys.executable)
    return app_path


def get_resources_file(fname):
    return os.path.join(get_application_path(), "resources", fname)


def load_ui_widget(uifilename, cls_to_register=None, parent=None):
    loader = QtUiTools.QUiLoader()
    loader.setLanguageChangeEnabled(True)
    loader.setTranslationEnabled(True)
    if cls_to_register:
        for cls in cls_to_register:
            loader.registerCustomWidget(cls)

    ui_fpath = os.path.join(get_application_path(), "forms", uifilename)
    if not os.path.isfile(ui_fpath):
        msg = '{} does not exist'.format(ui_fpath)
        raise ValueError(msg)

    print("Loading UI file: {0}".format(os.path.join(get_application_path(),
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

    frames_cnt = 0
    load_file = MDAnalysis.Universe(fname, format=format_data)
    if format_data == None:
        bonds = load_file.bonds.to_indices()
        no_bonds = len(load_file.bonds)
    return load_file, no_bonds
