# Scripts4Processing
Repository with useful scripts for SeisSol postprocessing.

In order to run the scripts, open Spyder (.py), Jupyter-Notebook (.ipynb) or MATLAB (.m) from the main directory in which the README file resides.

### How to clone this repository including its submodules?

Use
1. git clone --recurse-submodules https://github.com/fabian-kutschera/Scripts4Processing.git

If you already cloned this repository but forgot `--recurse-submodules` use:
1. git clone https://github.com/fabian-kutschera/Scripts4Processing.git
2. git submodule init
3. git submodule update

Find more information on how to use submodules in the [documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

### Overview of scripts


#### [concatenate_EnF_t.py](./concatenate_EnF_t.py)
This file is part of SeisSol. It is used to retrieve the Moment Magnitude Mw and the Scalar Seismic Moment M0, as well as to plot the Seismic Moment Rate and Frictional Energy as a function of time. The script requires an **EnF_t** output (e.g. *HFFtest-EnF_t-00000.dat*) of SeisSol.

#### [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py)
As an alternative to the [concatenate_EnF_t.py](./concatenate_EnF_t.py) one can use the computeMomentTensorSubSeisSol.py script. It requires the fault xdmf output (e.g. *HFFtest-fault.xdmf*). However, please provide a 3D velocity model in order to have an accurate estimate of the Moment Rate etc.. Remember to set the following flags:
- `--invertSls`: invert Slip in the strike direction due to SeisSol convention
- `--NZ 1` and `--NH 1` to plot the beachball representation (point source) later on
- `--asagiFile` or `--oneDvel` to include a 3D- or 1D-velocity model (3D is recommended)
- `--proj` to project back to WGS84 from given projection string (e.g. *--proj '+proj=utm +zone=27'* or *--proj '+init=EPSG:32627'*)
- `--outName` to append additional name to PointSourceFile***.h5 file (e.g. *--outName 'DR_complex_West'*)

The computeMomentTensorSubSeisSol.py generates a Point***.h5 output file which can be processed with several other script (see below).

#### [MomentRate_fromPointSource_h5.py](./MomentRate_fromPointSource_h5.py)
This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output. It will generate a **Moment Rate plot** similar to [concatenate_EnF_t.py](./concatenate_EnF_t.py).

TO DO!!!

#### [drawMapFromMomentTensorFile.py](./TeleseismicDataRelated/drawMapFromMomentTensorFile.py)
This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output. It will plot **the beachball on a map**. Consider the following flags for a decent map:

- `--proj`: **not** necessary if properly projected within [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py)
- TO DO!!!

#### [focal_plot.m](./focal_plot.m)
The `focal_plot.m` script requires to add [MATLAB code for moment tensor plotting](https://github.com/djpugh/MTplot) to the MATLAB path. The necessary files are already contained within this repository (MTplot) if cloned properly (see above).

This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output. It will generate a **detailed beachball plot** which should be consistent with the output from [drawMapFromMomentTensorFile.py](./TeleseismicDataRelated/drawMapFromMomentTensorFile.py).

TO DO!!!
