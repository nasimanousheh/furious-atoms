import os
import click
from setuptools import sandbox
import shutil
import PyInstaller.__main__


def create_installer(working_dir, use_debug=False):
    dist_dir = os.path.join(working_dir, 'installer_dist')
    build_dir = os.path.join(working_dir, 'installer_build')
    spec_dir = os.path.dirname(os.path.realpath(__file__))
    spec_file_rel = os.path.abspath(os.path.join(spec_dir,
                                                 'package_full.spec'))
    spec_file_debug = os.path.abspath(os.path.join(spec_dir,
                                                   'package_verbose.spec'))
    spec_file = spec_file_debug if use_debug else spec_file_rel
    PyInstaller.__main__.run(['--distpath', dist_dir,
                              '--workpath', build_dir,
                              spec_file])


def clean_installer(working_dir):
    dist_dir = os.path.join(working_dir, 'installer_dist')
    build_dir = os.path.join(working_dir, 'installer_build')
    for path in [dist_dir, build_dir]:
        if os.path.isdir(path):
            shutil.rmtree(path)


def clean_dist(working_dir):
    dist_dir = os.path.join(working_dir, 'dist')
    if os.path.isdir(dist_dir):
        shutil.rmtree(dist_dir)


def clean(working_dir):
    build_dir = os.path.join(working_dir, 'build')
    egg_dir = os.path.join(working_dir, 'furious_atoms.egg-info')
    for path in [egg_dir, build_dir]:
        if os.path.isdir(path):
            shutil.rmtree(path)


def install():
    sandbox.run_setup('setup.py', ['clean', 'install'])


def create_package():
    sandbox.run_setup('setup.py', ['clean', 'bdist_wheel'])


d_command = {"installer": create_installer,
             "clean-installer": clean_installer,
             "package": create_package,
             "clean-package": clean_dist,
             "clean": clean,
             "install": install
             }


@click.command()
@click.option('--verbose', '-v', is_flag=True, help="Print more output.")
@click.argument('name', type=click.Choice(d_command.keys()))
def main(name=None, verbose=False):
    # current_dir = os.getcwd()
    working_dir = os.path.dirname(os.path.realpath(__file__))
    working_dir = os.path.join(working_dir, os.pardir)
    name = name.lower()

    # os.chdir(working_dir)

    if name == 'installer':
        d_command.get(name)(working_dir, use_debug=verbose)
    elif name in ['clean-installer', 'clean-dist', 'clean']:
        d_command.get(name)(working_dir)
    else:
        d_command.get(name)()

    # os.chdir(current_dir)


if __name__ == "__main__":
    main()
