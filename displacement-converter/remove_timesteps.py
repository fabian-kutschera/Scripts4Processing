from netCDF4 import Dataset
import numpy as np
import argparse
from shutil import copyfile

parser = argparse.ArgumentParser(
    description="Remove timesteps from the end of a NetCDF file."
)
parser.add_argument("inputfile", help="netcdf input file")
parser.add_argument(
    "--timesteps",
    metavar="steps",
    type=int,
    default=2,
    help="number of timesteps to remove from the end of the file",
)
args = parser.parse_args()

inputfile = args.inputfile
outputfile = "modified_" + args.inputfile

print("Creating new NetCDF file with removed timesteps")

with Dataset(inputfile, mode="r") as src:
    with Dataset(outputfile, mode="w", format="NETCDF4") as dst:
        
        # Copy dimensions, making 'time' unlimited
        for name, dim in src.dimensions.items():
            dst.createDimension(name, None if name == "time" else len(dim))
        
        # Copy variables
        for name, var in src.variables.items():
            new_var = dst.createVariable(name, var.datatype, var.dimensions)
            new_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            
            # Copy data, removing the last 'steps_to_remove' timesteps for 'time' and 'z'
            if name == "time":
                dst.variables[name][:] = var[:-args.timesteps]
            elif name == "z":
                dst.variables[name][:] = var[:-args.timesteps, :, :]
            else:
                dst.variables[name][:] = var[:]

print("Modified NetCDF file saved as", outputfile)

