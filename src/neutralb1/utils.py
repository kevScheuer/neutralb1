"""A collection of some utility functions"""

import os
import pathlib
import re
import subprocess
import sys


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


def load_environment() -> None:
    """Load the shell environment variables from the .env file

    This function reads the `.env` file located in the `config` directory of the
    workspace directory, which is expected to be named `neutralb1`. It sets the
    environment variables defined in that file into the current Python environment.

    Raises:
        FileNotFoundError: If the `.env` file does not exist in the expected location.
    """

    workspace_dir = get_workspace_dir()

    if not os.path.exists(f"{workspace_dir}/config/.env"):
        raise FileNotFoundError(
            f"Environment file not found at {workspace_dir}/config/.env. "
            "Please execute `make update-env` to generate it."
        )

    env_vars = {}
    with open(f"{workspace_dir}/config/.env", "r") as env_file:
        for line in env_file:
            # Skip comments and empty lines
            if line.strip() and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                env_vars[key] = value

    os.environ.update(env_vars)
    sys.path.insert(0, workspace_dir)

    return


def get_workspace_dir() -> str:
    """Get the root directory of the workspace.

    Assumes the workspace directory is named "neutralb1" and the file calling this
    function is located somewhere within the workspace directory structure.

    Returns:
        str: The path of the workspace directory.
        Raises FileNotFoundError if the "neutralb1" directory is not found in the
    Raises:
        FileNotFoundError: If the "neutralb1" directory is not found in the path.
    """

    current_dir = pathlib.Path.cwd()
    workspace_dir = None
    while current_dir != current_dir.parent:
        if (current_dir / "neutralb1").is_dir():
            workspace_dir = str(current_dir / "neutralb1")

        current_dir = current_dir.parent
    if workspace_dir is None:
        raise FileNotFoundError("Could not find 'neutralb1' directory in the path.")

    return workspace_dir
