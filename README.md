# Script4Processing
Repository with useful scripts for SeisSol postprocessing.

#### How to clone this repository including it's submodules?

Use
1. git clone --recurse-submodules https://github.com/chaconinc/MainProject

If you already cloned this repository but forgot `--recurse-submodules` use:
1. git clone https://github.com/fabian-kutschera/Script4Processing.git
2. git submodule init
3. git submodule update

Find more information on how to use submodules in the [documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

#### Overview of scrips
The focal_plot.m(link) script requires to add [MATLAB code for moment tensor plotting](https://github.com/djpugh/MTplot) to the MATLAB path. The necessary files are already contained within this repository (foldername) if cloned properly (see above).
