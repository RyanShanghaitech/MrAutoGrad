from setuptools import setup, Extension
import numpy
import sys

modExt = Extension\
(
    "mrautograd.ext", 
    sources = 
    [
        './mrautograd_src/ext/mag/v3.cpp',
        './mrautograd_src/ext/mag/GradGen.cpp',
        './mrautograd_src/ext/mtg/mtg_functions.cpp',
        './mrautograd_src/ext/mtg/spline.cpp',
        './mrautograd_src/ext/main.cpp',
    ],
    libraries = [] if sys.platform=="win32" else ['jemalloc'],
    include_dirs = ["./mrautograd_src/ext/", numpy.get_include()],
    language = 'c++',
)

setup\
(
    name = 'mrautograd',
    install_requires = ["numpy", "matplotlib"],
    ext_modules = [modExt],
    packages = ["mrautograd", "mrautograd.ext", "mrautograd.ext.core", "mrautograd.ext.traj", "mrautograd.ext.mtg"],
    package_dir = {"mrautograd":"./mrautograd_src/", "mrautograd.ext":"./mrautograd_src/ext/", "mrautograd.ext.core":"./mrautograd_src/ext/mag/", "mrautograd.ext.traj":"./mrautograd_src/ext/traj/", "mrautograd.ext.mtg":"./mrautograd_src/ext/mtg/"},
    include_package_data = True
)
