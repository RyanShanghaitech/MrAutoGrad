# Magnetic Resonance Automatic Gradient Toolbox (MRAutoGrad, MAG)

## Usage
1. To install the pip package of the proposed algorithm (trajectory library and baseline algorithm is also included)
    ```
    $ pip install .
    ```
1. To reproduce the result in *A Graphical Method for Designing Time-Optimal Non-Cartesian Gradient Waveforms*, please run `benchmark.py`.
    ```
    $ python benchmark.py
    ```
1. To switch between different solver, edit this line in `benchmark.py`
    ```
    # select solver
    mag.setSolverMtg(0) # 0 for MAG solver, 1 for MTG solver
    ```
1. To switch between different trajectories, edit thes lines in `benchmark.py`, please leave one line uncommented and comment out other lines.
    ```
    def eval():
        global nAx
        # nAx = 2; return mag.getG_Spiral(bIs3D=0, **argCom) # 0.380s
        # nAx = 2; return mag.getG_VarDenSpiral(bIs3D=0, **argCom) # 0.499s
        nAx = 2; return mag.getG_Rosette(bIs3D=0, **argCom) # 16.39s (9 frames)
        # nAx = 2; return mag.getG_Rosette_Trad(**argCom, dOm1=10*pi, dOm2=8*pi, dTmax=1, dTacq=2e-03) # 16.39s (9 frames)
        # nAx = 3; return mag.getG_Shell3d(**argCom) # 183.7
        # nAx = 3; return mag.getG_Yarnball(**argCom) # 196.1
        # nAx = 3; return mag.getG_Seiffert(**argCom) # 232.9s
        # nAx = 3; return mag.getG_Cones(**argCom) # 149.9s
    ```
1. The generated gradient waveform and corresponding gspace curve, trajectory, slew rate and gradient amplitude are shown in the printed figure. The figure is also saved as `figure.pdf` in the root directory in case plotting is unavailable.
1. Number of iteration and execution time is shown in the command line interface, e.g:
    ```
    MAG Nit: 10122
    Elapsed time: 33.782 ms
    ```