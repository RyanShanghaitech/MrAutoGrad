from setuptools import setup, Extension
import numpy
from os.path import exists
from ctypes.util import find_library

useMtg = exists("./mrautograd_src/ext/mtg/")
useJemalloc = find_library("jemalloc") is not None

_sources = \
[
    './mrautograd_src/ext/utility/global.cpp',
    './mrautograd_src/ext/utility/v3.cpp',
    './mrautograd_src/ext/mag/GradGen.cpp',
    './mrautograd_src/ext/main.cpp',
    './mrautograd_src/ext/mtg/mtg_functions.cpp',
    './mrautograd_src/ext/mtg/spline.cpp'
]
if not useMtg:
    _sources.remove('./mrautograd_src/ext/mtg/mtg_functions.cpp')
    _sources.remove('./mrautograd_src/ext/mtg/spline.cpp')

modExt = Extension\
(
    "mrautograd.ext", 
    sources = _sources,
    libraries = ['jemalloc'] if useJemalloc else [],
    include_dirs = ["./mrautograd_src/ext/", numpy.get_include()],
    define_macros = [("USE_MTG", None)] if useMtg else None,
    language = 'c++'
)

_packages = \
[
    "mrautograd", 
    "mrautograd.ext", 
    "mrautograd.ext.traj",
    "mrautograd.ext.mag", 
    "mrautograd.ext.mtg",
    "mrautograd.ext.utility",
]
if not useMtg:
    _packages.remove("mrautograd.ext.mtg")

_package_dir = \
{
    "mrautograd":"./mrautograd_src/", 
    "mrautograd.ext":"./mrautograd_src/ext/", 
    "mrautograd.ext.traj":"./mrautograd_src/ext/traj/",
    "mrautograd.ext.mag":"./mrautograd_src/ext/mag/", 
    "mrautograd.ext.mtg":"./mrautograd_src/ext/mtg/",
    "mrautograd.ext.utility":"./mrautograd_src/ext/utility/"
}

setup\
(
    name = 'mrautograd',
    # install_requires = ["numpy", "matplotlib"], # pip will automatically upgrade numpy if it see this, which might corrupt the environment
    ext_modules = [modExt],
    packages = _packages,
    package_dir = _package_dir,
    include_package_data = True
)
