# -*- mode: python -*-

import os
import platform
from os.path import join as pjoin, dirname
import vtk
import numpy
import scipy
import MDAnalysis
from fury.data import DATA_DIR
import furiousatoms
import sys
sys.setrecursionlimit(5000)

block_cipher = None
VTK_PATH = dirname(vtk.__file__)
NUMPY_PATH = pjoin(dirname(numpy.__file__), ".libs")
SCIPY_PATH = pjoin(dirname(scipy.__file__), ".libs")
MDA_PATH = MDAnalysis.__path__ + [pjoin(dirname(MDAnalysis.__file__), 'lib'),
           pjoin(dirname(MDAnalysis.__file__), 'lib', 'format')]
FA_PATH = pjoin(dirname(furiousatoms.__file__), os.pardir)

FA_BIN = pjoin(FA_PATH, 'bin/furious-atoms')
if not os.path.isfile(FA_BIN):
    plt = platform.system()
    if plt == "Windows":
        FA_BIN = pjoin(os.path.dirname(sys.executable), 'Scripts/furious-atoms')
    else:
        FA_BIN = pjoin(os.path.dirname(sys.executable), 'furious-atoms')

added_files = [(pjoin(FA_PATH, 'furiousatoms/forms'), 'forms'),
               (pjoin(FA_PATH, 'furiousatoms/languages'), 'languages'),
               (pjoin(FA_PATH, 'furiousatoms/resources'), 'resources'),
               (pjoin(FA_PATH, 'furiousatoms/skybox0'), 'skybox0'),
               (DATA_DIR, 'fury/data/files'),
              ]

a = Analysis([FA_BIN],
             pathex=[VTK_PATH, NUMPY_PATH, SCIPY_PATH] + MDA_PATH,
             binaries=[],
             datas=added_files,
             hiddenimports=['cython', 'PySide2', 'PySide2.QtXml', 'shiboken2', 'scipy._lib.messagestream',
                            'vtkmodules', 'vtkmodules.all', 'vtkmodules.qt.QVTKRenderWindowInteractor',
                            'vtkmodules.util', 'vtkmodules.util.numpy_support', 'vtkmodules.numpy_interface',
                            'vtkmodules.numpy_interface.dataset_adapter', 'vtkmodules.util.colors',
                            'vtk', 'vtk.util', 'vtk.util.colors', 'MDAnalysis', 'MDAnalysis.lib.formats.cython_util',
                            'netCDF4_utils', 'netcdftime'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='furiousatoms',
          debug=False,
          strip=False,
          upx=True,
          console=True )

#######################################
# Code-sign the generated executable
# C:\Program Files (x86)\Windows Kits\10\bin\x64\signtool.exe
#import subprocess
#subprocess.call([
#   "SIGNTOOL.EXE",
#   "/F", "path-to-key.pfx",
#   "/P", "your-password",
#   "/T", "time-stamping url",
#   'name.exe',
#])
#######################################
