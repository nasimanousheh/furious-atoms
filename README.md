# Furious-Atoms


Small applications for atoms

## To install:

    pip install -e .

or

    python tools/make.py install

## To run it:

    furious-atoms

## To Create installer

    python tools/make.py installer
    python tools/make.py installer -v  # To get the verbose version

then you will the executable in `installer_dist` folder. To clean the build
you will have to run:

    python tools/make.py clean-installer

it will delete `installer_dist` and `installer_build`folders.

## To Create wheel/package

    python tools/make.py package

then, your wheel will be in `dist` folder. to clean it:

    python tools/make.py clean-package

## To clean

    python tools/make.py clean
