from setuptools import setup, Extension
import numpy
import sys, os

fMtgExist = os.path.exists('./mrautograd_src/ext/mtg')

_sources = \
[
    './mrautograd_src/ext/mag/v3.cpp',
    './mrautograd_src/ext/mag/GradGen.cpp',
    './mrautograd_src/ext/main.cpp',
]
if fMtgExist:
    _sources += \
    [
        './mrautograd_src/ext/mtg/mtg_functions.cpp',
        './mrautograd_src/ext/mtg/spline.cpp',
    ]

modExt = Extension\
(
    "mrautograd.ext", 
    sources = _sources,
    libraries = [] if sys.platform=="win32" else ['jemalloc'],
    include_dirs = ["./mrautograd_src/ext/", numpy.get_include()],
    language = 'c++',
    define_macros = [("MTG_EXIST",None)] if fMtgExist else []
)

_packages = \
[
    "mrautograd", 
    "mrautograd.ext", 
    "mrautograd.ext.core", 
    "mrautograd.ext.traj"
]
_package_dir = \
{
    "mrautograd":"./mrautograd_src/", 
    "mrautograd.ext":"./mrautograd_src/ext/", 
    "mrautograd.ext.core":"./mrautograd_src/ext/mag/", 
    "mrautograd.ext.traj":"./mrautograd_src/ext/traj/"
}
if fMtgExist:
    _packages += \
    [
        "mrautograd.ext.mtg"
    ]
    _package_dir["mrautograd.ext.mtg"] = "./mrautograd_src/ext/mtg/"

setup\
(
    name = 'mrautograd',
    install_requires = ["numpy", "matplotlib"],
    ext_modules = [modExt],
    packages = _packages,
    package_dir = _package_dir,
    include_package_data = True
)
