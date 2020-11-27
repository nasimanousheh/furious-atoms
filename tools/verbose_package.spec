# -*- mode: python -*-

import sys
import os
from os.path import join as pjoin

spec_root = os.path.abspath(os.path.join(SPECPATH, '..'))
print(spec_root)
sys.setrecursionlimit(5000)

block_cipher = None

added_files = [(pjoin(spec_root, 'furiousatoms/forms'), 'forms'),
               (pjoin(spec_root, 'furiousatoms/languages'), 'languages'),
               (pjoin(spec_root, 'furiousatoms/resources'), 'resources'),
              ]

a = Analysis([pjoin(spec_root, 'bin/furious-atoms')],
             pathex=[spec_root],
             binaries=[],
             datas=added_files,
             hiddenimports=['cython', 'sklearn', 'sklearn.neighbors.typedefs', 'PySide2', 'PySide2.QtXml',
                            'scipy._lib.messagestream', 'pkg_resources.py2_warn', 'MDAnalysis.lib.formats.cython_util'],
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
          exclude_binaries=True,
          name='FuriousAtoms',
          debug=False,
          strip=False,
          upx=False,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='FuriousAtoms')
