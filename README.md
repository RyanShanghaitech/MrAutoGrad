# Magnetic Resonance Automatic Gradient Toolbox (MRAutoGrad, MAG)

## Introduction
This toolbox is a pip package with C++ backend. The pip package can be called via Python interface to generate gradient waveforms for built-in and external trajectories. The C++ source code (in `mrautograd_src/ext/`) can be ported to other pulse sequence project like UIH's Adept project for gradient waveform calculation.

## Install
To install the pip package of the proposed algorithm (including trajectory library built on it):
```
$ pip install .
```
for linux user, `jemalloc` can be installed for improved performance. If `jemalloc` can not be installed, please remove this line in `setup.py`.
```
libraries = [] if sys.platform=="win32" else ['jemalloc'],
```

## Reproduce The Paper
1. To reproduce the result in *A Graphical Method for Designing Time-Optimal Non-Cartesian Gradient Waveforms*, please run `benchmark.py`.
    ```
    $ python benchmark.py
    ```
1. To switch between different solver, edit this line in `benchmark.py`.
    ```
    # select solver
    mag.setSolverMtg(0) # 0 for MAG solver, 1 for MTG solver
    ```
    Note that `MTG` solver might be unavailable in public repository due to copyright limit. The private repository which contiains `MTG` solver can be requested from the owner of this project via email.
1. To switch between different trajectories, edit these lines in `benchmark.py`, please leave one line uncommented and the other commented.
    ```
    def eval():
        global nAx
        # nAx = 2; return mag.getG_Spiral(bIs3D=0, **argCom)
        # nAx = 2; return mag.getG_VarDenSpiral(bIs3D=0, **argCom)
        nAx = 2; return mag.getG_Rosette(bIs3D=0, **argCom)
        # nAx = 2; return mag.getG_Rosette_Trad(**argCom, dOm1=10*pi, dOm2=8*pi, dTmax=1, dTacq=2e-03)
        # nAx = 3; return mag.getG_Shell3d(**argCom)
        # nAx = 3; return mag.getG_Yarnball(**argCom)
        # nAx = 3; return mag.getG_Seiffert(**argCom)
        # nAx = 3; return mag.getG_Cones(**argCom)
    ```
1. The generated gradient waveform and corresponding gspace curve, trajectory, slew rate and gradient amplitude are shown in the printed figure. The figure is also saved as `figure.pdf` in the root directory in case plotting is unavailable.
1. Number of iteration and execution time is shown in the command line interface, e.g:
    ```
    MAG Nit: 10122
    Elapsed time: 33.782 ms
    ```

## Examples
Examples for generating gradient waveforms for either built-in trajectory or external trajectory (expressed by trajectory function or trajectory samples) can be found in `example` folder.