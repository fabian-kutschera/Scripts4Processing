# Script4Processing
Repository with useful scripts for SeisSol postprocessing.

In order to run the scripts, open Spyder (.py), Jupyter-Notebook (.ipynb) or MATLAB (.m) from the main directory in which the README file resides.

### How to clone this repository including its submodules?

Use
1. git clone --recurse-submodules https://github.com/fabian-kutschera/Script4Processing.git

If you already cloned this repository but forgot `--recurse-submodules` use:
1. git clone https://github.com/fabian-kutschera/Script4Processing.git
2. git submodule init
3. git submodule update

Find more information on how to use submodules in the [documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

### Overview of scrips
#### [focal_plot.m](./focal_plot.m)
The `focal_plot.m` script requires to add [MATLAB code for moment tensor plotting](https://github.com/djpugh/MTplot) to the MATLAB path. The necessary files are already contained within this repository (MTplot) if cloned properly (see above).

#### [concatenate_EnF_t.py](./concatenate_EnF_t.py)
