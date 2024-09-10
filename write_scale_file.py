"""Write files to working directories to scale params in pseudo truth fits

The pseudo-truth cfg file requires a scale factor for the parameters. This script is a
quick method to write that scale factor to the working directories of the fits. This 
script could be included in the overall architecture of the analysis, but for now it is
standalone as pseudo-truth fits are not used often.
"""

import os

import pandas as pd

parent_dir = (
    "/home/kscheuer/volatile/TMPDIR/ampToolsFits/omegapi/"
    "allPeriods/PARA_0/ver03.1_mc/ver03/"
    "0m_1m_1p_iso/recoil-pi-mass_1.4/t_0.30-0.50/"
)
truth_df = pd.read_csv(
    "analysis/input-output-tests/0m_1m_1p/truth.csv", index_col="index"
)
scale_column = truth_df["par_scale"]

my_subdirs = []
for subdir, dirs, files in os.walk(parent_dir):
    if (
        "truth" in subdir
        and "log" not in subdir
        and "bootstrap" not in subdir
        and "rand" not in subdir
    ):
        my_subdirs.append(subdir)

my_subdirs.sort(key=lambda x: x.rsplit("-", 1)[1].split("/")[0])

for i, d in enumerate(my_subdirs):
    scale_factor = scale_column[i]
    with open(f"{d}/bestFitPars.txt", "w") as f:
        f.write(f"parameter par_scale {scale_factor} fixed")
