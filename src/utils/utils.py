"""A collection of some utility functions"""

import os
import pathlib
import re
import subprocess


def sort_input_files(input_files: list, position: int = -1) -> list:
    """Sort the input files based off the last number in the file name or path

    Args:
        input_files (list): input files to be sorted
        position (int, optional): Index position of the number to be sorted on in the
            full path. Defaults to -1, meaning the last number is used for sorting. Be
            careful using this, as it will assume all path names have the same amount of
            distinct numbers, and thus the same indices.

    Returns:
        list: sorted list of files
    """

    def extract_last_number(full_path: str) -> float:
        numbers = re.findall(r"(?:\d*\.*\d+)", full_path)
        return float(numbers[position]) if numbers else float("inf")

    return sorted(input_files, key=extract_last_number)


def load_shell_environment() -> None:
    """Load the shell environment variables from the setup script.

    This function executes the `setup_gluex.sh` script in a bash shell and captures the
    environment variables it sets. It then updates the current Python environment with
    these variables. Particularly useful for jupyter notebook sessions.
    """

    # Load the environment variables from the shell setup script
    root_dir = pathlib.Path(__file__).parent.parent.parent
    command = f"bash -l -c 'source {root_dir}/setup_gluex.sh && env'"
    proc = subprocess.Popen(
        command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
    )
    output, _ = proc.communicate()
    # Parse the environment variables
    env_vars = {}
    for line in output.decode().splitlines():
        # the output contains a bunch of BASH_FUNCS that will ruin the environment
        # variables. this avoids those issues
        if (
            len(line.split("=", 1)) != 2
            or line.startswith("BASH_FUNC")
            or line.startswith(" ")
            or line.startswith("\t")
        ):
            continue
        key, value = line.split("=", 1)
        env_vars[key] = value
    os.environ.update(env_vars)

    return
