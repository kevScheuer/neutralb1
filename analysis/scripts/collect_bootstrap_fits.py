""" Collect all bootstrap fits in a t-bin into a csv, labelled by their mass bin

This script will iterate over all the mass subdirectories of a t-bin, and collect their
bootstrap fits. Every row of the resulting csv is labelled with the mass range it was 
done in, and the bin index it corresponds to. The mass ranges are also sorted so the 
csv is easier to read directly
"""

import os

import pandas as pd


def main(parent_dir: str):
    # first check if ROOT is loaded into the shell environment
    if not os.environ["ROOTSYS"]:
        raise EnvironmentError("ROOTSYS path is not loaded\n")

    # get all the mass bins and store their ranges for sorting and labelling
    bin_ranges = []
    for subdir, dirs, files in os.walk(parent_dir):
        # ensure subdir always contains the [rand, bootstrap] subdirectories
        if "rand" not in dirs or "bootstrap" not in dirs:
            continue

        # hardcoded to grab "binName_#-#" style of subdirectory
        t_range = subdir.split("/")[-2]
        bin_ranges.append(subdir.split("/")[-1].split("_")[-1])

    bin_ranges = sorted(bin_ranges)

    # now we loop over the sorted list
    df_list = []
    for i, range in enumerate(bin_ranges):
        bootstrap_dir = f"{parent_dir}mass_{range}/bootstrap/"
        os.system(
            (
                "root -l -b -q"
                " 'analysis/scripts/fitsToCsv.C(\""
                f"-d {bootstrap_dir}"
                f" -o {range}.csv"
                "\")'"
            )
        )

        # add the bin index and range these fits are associated with
        df = pd.read_csv(f"{range}.csv", index_col="index")
        df["bin_range"] = range
        df["bin"] = i

        df_list.append(df)
        os.remove(f"{range}.csv")  # remove file so working dir isn't cluttered

    pd.concat(df_list).to_csv(f"bootstrap-fits_{t_range}")

    pass


if __name__ == "__main__":
    # TODO: build out argparse for these variables
    parent_dir = (
        "/lustre19/expphy/volatile/halld/home/kscheuer/ampToolsFits/omegapi"
        "/allPeriods/PARA_0/ver03.1_mc/ver03/0m_1m_1p_iso"
        "/recoil-pi-mass_1.4/t_0.30-0.50/"
    )
    main(parent_dir)
