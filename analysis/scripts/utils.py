"""A collection of some utility functions used in the analysis scripts."""

import re


def sort_input_files(input_files: list) -> list:
    """Sort the input files based off the last number in the file name or path"""

    def extract_last_number(full_path: str) -> float:
        numbers = re.findall(r"(?:\d*\.*\d+)", full_path)
        return float(numbers[-1]) if numbers else float("inf")

    return sorted(input_files, key=extract_last_number)
