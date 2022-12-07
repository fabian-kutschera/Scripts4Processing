import h5py
from tqdm import tqdm
import re
import os
import shutil
import argparse

dic = {"u1": "U", "u2": "V", "u3": "W"}

parser = argparse.ArgumentParser(description="update surface file u1u2u3 to UVW")
parser.add_argument("xdmf_filename", help="filename")
args = parser.parse_args()

fn = args.xdmf_filename
fnew = fn.split("-surface.xdmf")[0] + "UVW-surface.xdmf"

print(f"modifying variable names in {fnew} ...")
with open(fn, "r") as sources:
    lines = sources.readlines()
assert fn != fnew

with open(fnew, "w") as sources:
    for line in lines:
        line_mod = line
        line_mod = re.sub(r"mesh0/u1", "mesh0/U", line_mod)
        line_mod = re.sub(r"mesh0/u2", "mesh0/V", line_mod)
        line_mod = re.sub(r"mesh0/u3", "mesh0/W", line_mod)
        line_mod = re.sub(r"Name=\"u1\"", 'Name="U"', line_mod)
        line_mod = re.sub(r"Name=\"u2\"", 'Name="V"', line_mod)
        line_mod = re.sub(r"Name=\"u3\"", 'Name="W"', line_mod)
        line_mod = re.sub(r"-surface_cell.h5", "UVW-surface_cell.h5", line_mod)
        sources.write(line_mod)

fnh5 = os.path.splitext(fn)[0] + "_cell.h5"
fnewh5 = os.path.splitext(fnew)[0] + "_cell.h5"
print(f"copying {fnh5} to {fnewh5}")
if not os.path.exists(fnh5):
    raise FileNotFoundError(fnh5)

shutil.copyfile(fnh5, fnewh5)

print(f"modifying variable names in {fnewh5}")
with h5py.File(fnewh5, "r+") as f:
    # Iterate over each group
    for top_key, group in f.items():
        # Rename all datasets in the group (pad with zeros)
        for key in tqdm(group.keys()):
            if key in dic:
                new_key = dic[key]
                group.move(key, new_key)
