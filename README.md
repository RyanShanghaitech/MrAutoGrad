# Magnetic Resonance Automatic Gradient Toolbox (MRAutoGrad, MAG)

## Introduction
This toolbox is a pip package with C++ backend. The pip package can be called via Python interface to generate **non-Cartesian** gradient waveforms for built-in and external trajectories. The C++ source code (in `mrautograd_src/ext/`) can be ported to other pulse sequence project like UIH's Adept project for gradient waveform calculation.

## Install
**Optionally**, to create a new conda environment (in case the dependencies in this package break your current environment), please run:
```bash
$ conda create -n magtest -y
$ conda activate magtest
```

To install the pip package of the proposed algorithm (including trajectory library built on it), and also the dependencies:
```bash
$ bash install.bash
```

It's the best practice to use my script `install.bash` for installation. You can also install via `pip install .` but remember to delete `*.egg-info` or pip will run into bug when uninstalling this package in current folder (see comments in `install.bash`).  

For linux user, optionally, `jemalloc` can be installed for improved performance.  

## Examples & Usages
Examples for generating gradient waveforms for either built-in trajectory (trajectory library) or external trajectory (expressed by trajectory function or trajectory samples) can be found in the `example` folder.

## Citation
If this project helps you, please cite [our paper](https://arxiv.org/abs/2507.21625):

[1] R. Luo, H. Huang, Q. Miao, J. Xu, P. Hu, and H. Qi, “Real-Time Gradient Waveform Design for Arbitrary k-Space Trajectories,” Sep 9, 2025, arXiv preprint arXiv:2507.21625.


