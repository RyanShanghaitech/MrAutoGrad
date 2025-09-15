# Magnetic Resonance Automatic Gradient Toolbox (MRAutoGrad, MAG)

## Introduction
This toolbox is a pip package with C++ backend. The pip package can be called via Python interface to generate **non-Cartesian** gradient waveforms for built-in and external trajectories. The C++ source code (in `mrautograd_src/ext/`) can be ported to other pulse sequence project like UIH's Adept project for gradient waveform calculation.

## Install
To install the pip package of the proposed algorithm (including trajectory library built on it):
```
$ bash install.bash
```
It's the best practice to use my script `install.bash` for installation. You can also install via `pip install .` but remember to delete `*.egg-info` or pip will run into bug when uninstalling this package in current folder (`*.egg-info` will mislead the pip about the package location).  

For linux user, optionally, `jemalloc` can be installed for improved performance.  

`numpy` is needed as a dependency.  

## Examples & Usages
Examples for generating gradient waveforms for either built-in trajectory or external trajectory (expressed by trajectory function or trajectory samples) can be found in `example` folder.

## Citation
If this project helps you, please cite [our paper](https://arxiv.org/abs/2507.21625):

[1] R. Luo, H. Huang, Q. Miao, J. Xu, P. Hu, and H. Qi, “Real-Time Gradient Waveform Design for Arbitrary k-Space Trajectories,” Sep 9, 2025, arXiv preprint arXiv:2507.21625.


