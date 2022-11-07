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

Run from within Spyder via `run concatenate_EnF_t.py /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/complex_fault_geometry/Complex_Middle_M7.07/HFFtest --plot`

#### [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py)
As an alternative to the [concatenate_EnF_t.py](./concatenate_EnF_t.py) one can use the computeMomentTensorSubSeisSol.py script. It requires the fault xdmf output (e.g. *HFFtest-fault.xdmf*). However, please provide a 3D velocity model in order to have an accurate estimate of the Moment Rate etc.. Remember to set the following flags:
- `--invertSls`: invert Slip in the strike direction due to SeisSol convention
- `--NZ 1` and `--NH 1` to plot the beachball representation (point source) later on
- `--asagiFile` or `--oneDvel` to include a 3D- or 1D-velocity model (3D is recommended)
- `--proj` to project back to WGS84 from given projection string (e.g. *--proj '+proj=utm +zone=27'* or *--proj '+init=EPSG:32627'*)
- `--outName` to append additional name to PointSourceFile***.h5 file (e.g. *--outName 'DR_complex_West'*)

The computeMomentTensorSubSeisSol.py generates a Point***.h5 output file which can be processed with several other scripts (see below).

In case you encounter any issues with this script, please consult [TeleseismicDataRelated](https://gitlab.lrz.de/thomas.ulrich/TuSeisSolScripts/-/tree/master/TeleseismicDataRelated) and use instead [compute_multi_cmt.py](https://gitlab.lrz.de/thomas.ulrich/TuSeisSolScripts/-/blob/master/TeleseismicDataRelated/compute_multi_cmt.py).

#### [MomentRate_fromPointSource_h5.py](./MomentRate_fromPointSource_h5.py)
This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output; adjust the input *fname* accordingly. It will generate a **Moment Rate Plot** similar to [concatenate_EnF_t.py](./concatenate_EnF_t.py).

Note: In order to obtain an accurate Moment Rate use the `--asagiFile` flag within [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py).

#### [drawMapFromMomentTensorFile.py](./TeleseismicDataRelated/drawMapFromMomentTensorFile.py)
This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output. It will plot **the beachball on a map**. Run the script from within Spyder e.g. via 

`run drawMapFromMomentTensorFile.py ../output/PointSourceFile_1_1_Simple_Middle_M7.333.h5 --MapBoundaries -20 -16 65 67 --outName 'Simple_Middle_M7.333' --Title 'Simple Middle M7.33'`

The flag `--proj` is **not** necessary if properly projected within [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py).

In case you encounter any issues with this script, please consult [TeleseismicDataRelated](https://gitlab.lrz.de/thomas.ulrich/TuSeisSolScripts/-/tree/master/TeleseismicDataRelated) and use instead [plot_map_cmt.py](https://gitlab.lrz.de/thomas.ulrich/TuSeisSolScripts/-/blob/master/TeleseismicDataRelated/plot_map_cmt.py).


#### [focal_plot.m](./focal_plot.m)
The `focal_plot.m` script requires to add [MATLAB code for moment tensor plotting](https://github.com/djpugh/MTplot) to the MATLAB path. The necessary files are already contained within this repository (MTplot) if cloned properly (see above).

This script can be used after [computeMomentTensorSubSeisSol.py](./TeleseismicDataRelated/computeMomentTensorSubSeisSol.py) since it requires the *PointSourceFile.h5* output. It will generate a **detailed beachball plot** which should be consistent with the output from [drawMapFromMomentTensorFile.py](./TeleseismicDataRelated/drawMapFromMomentTensorFile.py).

The [focal_plot.m](./focal_plot.m) finds all *PointSource.h5* files within the given directory and saves the corresponding beachballs seperately as .jpg figures.

#### [ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py](./ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py)

A succesful run of the script requires several installations, like the pythonXDMFReader as well as gmpe-smtk + latest openquake engine.
Procede as follows:

Download gmpe-smtk and pythonXdmfReader module
1. git clone https://github.com/GEMScienceTools/gmpe-smtk 
2. git clone https://gitlab.lrz.de/thomas.ulrich/pythonXdmfReader.git or pip install seissolxdmf via https://pypi.org/project/seissolxdmf/

Create links to smtk and pythonXdmfReader folders
- ln -s $LOCAL_DIR/gmpe-smtk/smtk
- ln -s $LOCAL_DIR/pythonXdmfReader

Run python code (on single processor):

- `python ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py --noMPI Sulawesi-surface.xdmf`   
or a faster run (parallel with MPI) , here assuming 4 cores processors 
- `python ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py --MP 4 Sulawesi-surface.xdmf`  
- e.g. `python ../../../Scripts/Script4Processing/ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py --MP 4 HFFtest-surface.xdmf`

Additional packages might need to be installed; as for example [Mpi4Py](https://anaconda.org/conda-forge/mpi4py) to run the code in parallel with MPI.

#### [CreateVtkCoastLineFromGmt.py](./CreateVtkCoastLineFromGmt.py)

Generates a .vtk file of the coastline using GMT. The CoastLine.vtk file can be opened using ParaView. The original script can be viewed in the official [SeisSol repository](https://github.com/SeisSol/SeisSol/blob/master/postprocessing/visualization/tools/CreateVtkCoastLineFromGmt.py).

Run the script from within Spyder (on heisenbug) e.g. via 

`run CreateVtkCoastLineFromGmt.py --lon -20.5 -15.5 --lat 65 67 --proj '+proj=utm +zone=27' --resolution "h"`

#### [compute_diff_seissol_data.py](./compute_diff_seissol_data.py)

Make difference between 2 (paraview) output files: f2-f1. The output must be from the same mesh, but the partionning may differ.

Run the script e.g. via

`python /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/SeisSol/postprocessing/science/compute_diff_seissol_data.py output_o4/HFFZ_fullycp_o4-surface.xdmf output_o6/HFFZ_fullycp_o6-surface.xdmf --idt 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 --Data u3 v3`

#### [extractDataFromUnstructuredOutput.py](extractDataFromUnstructuredOutput.py)

This script can help to resample SeisSol .xdmf output files and/ or extract timesteps/ variables of interest. The script depends on [recreateXdmf.py](recreateXdmf.py), obtained from [SeisSol](https://github.com/SeisSol/SeisSol/tree/master/postprocessing/visualization/tools).

Run the script from the terminal e.g. via

`python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/extractDataFromUnstructuredOutput.py HFFZ_fullycp_o4_300s-surface.xdmf --idt $(seq 0 720) --Data u3 locationFlag`

#### [vizualizeBoundaryConditions.py](vizualizeBoundaryConditions.py)

Double check that the boundary conditions on the mesh are correct with e.g.:

`python ../../Scripts4Processing/vizualizeBoundaryConditions.py Samos_WL.xdmf 0`
