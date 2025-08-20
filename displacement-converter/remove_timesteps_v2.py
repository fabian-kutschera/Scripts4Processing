#!/usr/bin/env python3
from netCDF4 import Dataset
import argparse
import os
import sys

def main():
    p = argparse.ArgumentParser(
        description="Trim a NetCDF file by keeping the first M time steps (or removing the last N)."
    )
    p.add_argument("inputfile", help="NetCDF input file")
    p.add_argument("-o", "--output", default=None, help="Output file (default: trimmed_<inputfile>)")
    p.add_argument("--keep", type=int, default=None,
                   help="Keep the first M time steps (preferred option)")
    p.add_argument("--remove", "--timesteps", dest="remove", type=int, default=None,
                   help="Remove the last N time steps (alias: --timesteps)")
    p.add_argument("--time-dim", default="time",
                   help="Name of the time dimension (default: 'time')")
    args = p.parse_args()

    infile = args.inputfile
    outfile = args.output or f"modified_{os.path.basename(infile)}"

    with Dataset(infile, "r") as src, Dataset(outfile, "w", format=src.file_format) as dst:
        # Copy global attributes
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})

        # Ensure time dimension exists
        if args.time_dim not in src.dimensions:
            sys.exit(f"ERROR: Dimension '{args.time_dim}' not found in {infile}")

        nsteps = len(src.dimensions[args.time_dim])

        # Determine how many to keep
        if args.keep is not None:
            keep = args.keep
        else:
            remove = args.remove if args.remove is not None else 1  # backward-compatible default
            keep = nsteps - remove

        if keep < 0:
            sys.exit(f"ERROR: Requested keep={keep} < 0 (nsteps={nsteps}).")
        if keep > nsteps:
            keep = nsteps  # cap to available steps

        # Create dimensions (force time to be unlimited)
        for name, dim in src.dimensions.items():
            if name == args.time_dim:
                dst.createDimension(name, None)  # unlimited
            else:
                dst.createDimension(name, len(dim))

        # Create and copy variables
        for name, var in src.variables.items():
            # Preserve compression/chunking where possible
            kw = {}
            try:
                f = getattr(var, "filters", None)
                if callable(f):
                    filt = var.filters()
                    if filt.get("zlib"):
                        kw.update(zlib=True, complevel=filt.get("complevel", 4))
                chunks = var.chunking()
                if isinstance(chunks, tuple):
                    kw["chunksizes"] = chunks
            except Exception:
                pass

            # Handle _FillValue
            fill_val = getattr(var, "_FillValue", None)
            if fill_val is not None:
                kw["fill_value"] = fill_val

            new_var = dst.createVariable(name, var.datatype, var.dimensions, **kw)

            # Copy variable attributes (but avoid duplicating _FillValue)
            new_var.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})

            # Data copy
            if args.time_dim in var.dimensions:
                if keep == 0:
                    # Nothing to write; unlimited dim length will remain 0 until something is written
                    continue
                t_axis = var.dimensions.index(args.time_dim)

                # Build a slice that selects the first `keep` steps along the time axis
                slicer = [slice(None)] * var.ndim
                slicer[t_axis] = slice(0, keep)

                # IMPORTANT: also slice on the LHS to avoid broadcasting mismatch
                new_var[tuple(slicer)] = var[tuple(slicer)]
            else:
                new_var[:] = var[:]

    print(f"Wrote {outfile} with {keep}/{nsteps} time steps kept.")

if __name__ == "__main__":
    main()

