"""A collection of methods for preprocessing PWA fit result DataFrames."""

import cmath
import pathlib
import warnings
from typing import Optional

import numpy as np
import pandas as pd

import neutralb1.analysis.physics as physics
import neutralb1.utils as utils


def standardize_fit_types(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize the types of the columns in a fit result DataFrame.

    All fit result dataframes have common columns that need not be large float values.
    This function converts those columns to smaller types to save memory.

    Args:
        df (pd.DataFrame): The DataFrame to standardize.
    Returns:
        pd.DataFrame: The DataFrame with standardized types.
    """
    info_columns_to_types = {
        "eMatrixStatus": "int8",
        "lastMinuitCommandStatus": "int8",
        "file": "category",
    }
    zero_columns = [col for col in df.columns if (df[col] == 0).all()]
    phase_columns = utils.get_phase_differences(df)

    # warn the user if phase columns are not yet wrapped
    if any(df[col].max() > 180.0 or df[col].min() < -180.0 for col in phase_columns):
        warnings.warn(
            "The following columns contain phases that are not wrapped to (-180, 180]: "
            f"{', '.join(phase_columns)}. Consider wrapping them before proceeding.",
            UserWarning,
        )

    return (
        df.astype(info_columns_to_types)
        .astype({c: "int8" for c in zero_columns} if zero_columns else {})
        .astype({c: "float32" for c in phase_columns} if phase_columns else {})
    )


def standardize_data_types(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize the types of the columns in a data file DataFrame.

    A 'data' DataFrame contains the values of the ROOT file that is being fit to.
    This function converts the columns to smaller types to save memory.

    Args:
        df (pd.DataFrame): The DataFrame to standardize.
    Returns:
        pd.DataFrame: The DataFrame with standardized types.
    """
    # the low/high/center columns are for TEM bins, which do not require precision or
    # bounds larger than float16
    f16_columns = [c for c in df.columns if "low" in c or "high" in c or "center" in c]
    # float32 should be the maximum for the average, RMS, and event count columns
    f32_columns = [c for c in df.columns if "avg" in c or "rms" in c or "events" in c]

    columns_to_types = {
        "file": "category",
        **{c: "float16" for c in f16_columns if f16_columns},
        **{c: "float32" for c in f32_columns if f32_columns},
    }

    return df.astype(columns_to_types)


def standardize_moment_types(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize the types of the columns in a moment DataFrame.

    A 'moment' DataFrame contains the values of the fitted moments, or moment projected
    from the fit results. This function converts the columns to smaller types to save
    memory.

    Args:
        df (pd.DataFrame): The DataFrame to standardize.
    Returns:
        pd.DataFrame: The DataFrame with standardized types.
    """
    moment_columns = [c for c in df.columns if c.startswith("H")]

    columns_to_types = {"file": "category", **{c: "float32" for c in moment_columns}}

    return df.astype(columns_to_types)


def find_null_columns(df: pd.DataFrame) -> list:
    """Find columns in a DataFrame that contain null values.

    Args:
        df (pd.DataFrame): The DataFrame to check for null values.
    Returns:
        list: A list of column names that contain null values.
    """
    return df.columns[df.isnull().any()].tolist()


def link_dataframes(
    fit_df: pd.DataFrame, target_df: pd.DataFrame, max_depth: int = 2
) -> pd.DataFrame:
    """Add a 'fit_index' to a target DataFrame to track which fit it belongs to.

    Randomized, bootstrap, truth, and data DataFrames have "file" columns that denote
    the path to the file each result row was derived from. The files are first
    expected to be subdirectories of the fit DataFrame's "file" column, but if not, then
    the function will search for the closest parent directory that matches, up to a
    maximum depth of `max_depth`. If found, the function will add a 'fit_index' column
    to the target DataFrame, which indicates the index of the fit DataFrame that
    corresponds to each file in the target DataFrame.

    Args:
        fit_df (pd.DataFrame): The fit DataFrame containing file paths.
        target_df (pd.DataFrame): The target DataFrame to link to. It's file columns
            must be subdirectories of the fit DataFrame's file column.
        max_depth (int, optional): The maximum depth to search for a matching fit file.
            Defaults to 2, meaning it will look for a fit file that is, or has a parent
            directory that is, at most 2 levels up from the target file's directory.
    Returns:
        pd.DataFrame: The target_df with an added 'fit_index' column that indicates
            the index of the fit DataFrame that corresponds to each file in the target
            DataFrame.

    Raises:
        KeyError: If the 'file' column is not present in either DataFrame.
        ValueError: If no matching fit is found for a file in the target DataFrame.
    """
    if "file" not in fit_df.columns:
        raise KeyError("The fit DataFrame must contain a 'file' column.")
    if "file" not in target_df.columns:
        raise KeyError("The target DataFrame must contain a 'file' column.")

    # Pre-resolve and cache all fit file paths and their parent hierarchies
    fit_path_cache = {}
    for idx, file_str in enumerate(fit_df["file"].astype(str)):
        resolved = pathlib.Path(file_str).resolve()
        # Store ancestors up to max_depth
        ancestors = []
        current = resolved.parent
        for _ in range(max_depth + 1):
            ancestors.append(current)
            if current == current.parent:
                break
            current = current.parent
        fit_path_cache[idx] = (resolved, set(ancestors))

    # Build a lookup dict: ancestor_path -> list of (fit_index, depth)
    ancestor_to_fits = {}
    for fit_idx, (fit_path, ancestors) in fit_path_cache.items():
        for depth, ancestor in enumerate(ancestors):
            if depth <= max_depth:
                if ancestor not in ancestor_to_fits:
                    ancestor_to_fits[ancestor] = []
                ancestor_to_fits[ancestor].append((fit_idx, depth))

    def find_fit_index(target_file: str) -> int:
        """If target_path shares a common ancestor with a fit file, return its index."""
        target_path = pathlib.Path(target_file).resolve()

        # Get target ancestors up to max_depth
        target_ancestors = []
        current = target_path.parent
        for depth in range(max_depth + 1):
            target_ancestors.append((current, depth))
            if current == current.parent:
                break
            current = current.parent

        # Find first matching ancestor
        for target_ancestor, target_depth in target_ancestors:
            if target_ancestor in ancestor_to_fits:
                # Return the first fit that matches within depth constraints
                for fit_idx, fit_depth in ancestor_to_fits[target_ancestor]:
                    if fit_depth <= max_depth and target_depth <= max_depth:
                        return fit_idx

        raise FileNotFoundError(f"No matching fit found for {target_file}")

    # Get unique file paths from target DataFrame to avoid redundant lookups
    unique_files = target_df["file"].astype(str).unique()

    # Build mapping for unique files only
    file_to_fit_index = {file: find_fit_index(file) for file in unique_files}

    # Map all rows at once using the pre-computed mapping
    return target_df.assign(
        fit_index=target_df["file"].astype(str).map(file_to_fit_index).astype("uint16")
    )


def add_missing_columns(
    reference_df: pd.DataFrame, target_df: pd.DataFrame
) -> pd.DataFrame:
    """Add column of zeros to the target if it is missing columns from the reference.

    Args:
        reference_df (pd.DataFrame): The reference DataFrame to check against.
        target_df (pd.DataFrame): The target DataFrame to modify.

    Returns:
        pd.DataFrame: The modified target DataFrame with missing columns added.
    """
    missing_cols = set(reference_df.columns) - set(target_df.columns)
    if not missing_cols:
        return target_df

    # Create DataFrame with missing columns all at once to avoid fragmentation
    missing_data = pd.DataFrame(
        {col: pd.Series(0, index=target_df.index, dtype="int8") for col in missing_cols}
    )

    return pd.concat([target_df, missing_data], axis=1)


def wrap_radians(num: float) -> float:
    """Wrap phase differences to be from (-pi, pi] & convert from radians to degrees

    Args:
        num (float): The phase to be wrapped
    Returns:
        float: The wrapped phase in degrees
    """
    return np.rad2deg(np.angle(np.exp(1j * num)))


def wrap_phases(df: pd.DataFrame, phase_columns: Optional[list] = None) -> pd.DataFrame:
    """Wrap phase columns in a DataFrame to be from (-pi, pi].

    Args:
        df (pd.DataFrame): The DataFrame containing phase columns to be wrapped.
            Expects the columns to be in radians.
        phase_columns (Optional[list], optional): A list of phase column names to wrap.
            Defaults to None, which will find the phase columns using
            :func:`neutralb1.utils.get_phase_differences`.

    Returns:
        pd.DataFrame: The DataFrame with wrapped phase columns.
    """
    if phase_columns is None:
        phase_columns = list(set(utils.get_phase_differences(df)))
        phase_columns += [f"{col}_err" for col in phase_columns]
    if "dphase" in df.columns:
        phase_columns.append("dphase")
        phase_columns.append("dphase_err")

    return df.assign(**{col: df[col].map(wrap_radians) for col in phase_columns})


def align_phase_difference_names(
    reference_df: pd.DataFrame,
    target_df: pd.DataFrame,
    phase_columns: Optional[list] = None,
) -> pd.DataFrame:
    """Rename the phase differences of a target to match those in the source DataFrame.

    Phase differences are named in the format 'eJPmL1_eJPmL2', but can have reverse
    ordering between two DataFrames. This function renames the columns in the target
    DataFrame to match the names in the source DataFrame, ensuring that the phase
    differences are consistent across both DataFrames.

    Args:
        reference_df (pd.DataFrame): The DataFrame whose phase columns will be used as
            the point of reference for naming.
        target_df (pd.DataFrame): The DataFrame whose phase columns will be renamed to
            match the reference DataFrame.
        phase_columns (Optional[list], optional): A list of phase column names to align.
            Defaults to None, which will find the phase columns using
            :func:`neutralb1.utils.get_phase_differences`.

    Returns:
        pd.DataFrame: The DataFrame with aligned phase columns.
    """
    if phase_columns is None:
        phase_columns = list(set(utils.get_phase_differences(reference_df)))

    reversed_phase_columns = ["_".join(pd.split("_")[::-1]) for pd in phase_columns]

    return target_df.rename(
        columns={
            old: new
            for old, new in zip(reversed_phase_columns, phase_columns)
            if old in target_df.columns and new in reference_df.columns
        }
    )


def restore_breit_wigner_phases(
    df: pd.DataFrame, mass_bins: pd.Series, phase_columns: Optional[list] = None
) -> pd.DataFrame:
    """Adjust the phases of a dataframe with breit-wigner amplitudes to correct values.

    A DataFrame produced from a fit result that contains breit-wigner amplitudes will
    not have the correct phase differences, due to how AmpTools calculates the phase
    difference. This function is used to restore the phases of such DataFrames to their
    correct values.

    Note:
        The error for the phase differences will no longer be meaningful, but this
        application is normally done for truth DataFrames, which do not have errors.

    Args:
        df (pd.DataFrame): The DataFrame containing truth phase columns to be reset.
        mass_bins (pd.Series): The center of the mass bins for each fit result row.
        phase_columns (Optional[list], optional): A list of phase column names to reset.
            Defaults to None, and will find the phase columns using
            :func:`neutralb1.utils.get_phase_differences`.

    Returns:
        pd.DataFrame: The DataFrame with reset truth phase columns.
    """

    def new_phase_dif(
        mass: float,
        amp1_re: float,
        amp1_im: float,
        bw_mass1: float,
        bw_width1: float,
        bw_l1: int,
        amp2_re: float,
        amp2_im: float,
        bw_mass2: float,
        bw_width2: float,
        bw_l2: int,
    ):
        complex_val1 = complex(amp1_re, amp1_im)
        complex_val2 = complex(amp2_re, amp2_im)
        bw1 = physics.breit_wigner(mass, bw_mass1, bw_width1, bw_l1)
        bw2 = physics.breit_wigner(mass, bw_mass2, bw_width2, bw_l2)
        phase1 = cmath.phase(bw1 * complex_val1)
        phase2 = cmath.phase(bw2 * complex_val2)

        return phase1 - phase2

    if phase_columns is None:
        phase_columns = list(set(utils.get_phase_differences(df)))

    df = df.copy()
    for pd in phase_columns:
        amp1, amp2 = pd.split("_")

        jp1 = amp1[1:3]
        l1 = utils.char_to_int(amp1[-1])
        jp2 = amp2[1:3]
        l2 = utils.char_to_int(amp2[-1])

        df[pd] = np.vectorize(new_phase_dif)(
            mass_bins,
            df[f"{amp1}_re"],
            df[f"{amp1}_im"],
            df[f"{jp1}_mass"],
            df[f"{jp1}_width"],
            l1,
            df[f"{amp2}_re"],
            df[f"{amp2}_im"],
            df[f"{jp2}_mass"],
            df[f"{jp2}_width"],
            l2,
        )

    return df


def filter_projected_moments(df: pd.DataFrame) -> pd.DataFrame:
    """Remove projected moment columns that are expected to be zero.

    When projecting the moments from the fit results, all moments are split
    into their real and imaginary parts. The H0 and H1 moments should be purely real,
    and the H2 moment should be purely imaginary. This function checks that this is
    the case, and removes H0 and H1 imaginary parts, and H2 real parts.

    Args:
        df (pd.DataFrame): The DataFrame containing projected moment columns.
    Returns:
        pd.DataFrame: The DataFrame with the specified projected moment columns removed.
    """
    cols_to_remove = []
    for col in df.columns:
        if (col.startswith("H0") or col.startswith("H1")) and col.endswith("_imag"):
            cols_to_remove.append(col)
        elif col.startswith("H2") and col.endswith("_real"):
            cols_to_remove.append(col)

    for col in cols_to_remove:
        if not np.allclose(df[col], 0, atol=1e-3):
            raise ValueError(
                f"Projected moment column {col} is expected to be 0 but contains"
                " non-zero values."
            )

    return df.drop(columns=cols_to_remove)


def remove_real_imag_suffixes(df: pd.DataFrame) -> pd.DataFrame:
    """Remove the '_real' and '_imag' suffixes from projected moment columns.

    .. note::
        This function should be called after :func:`filter_projected_moments`, which
        removes the projected moment columns that are expected to be zero.

    Args:
        df (pd.DataFrame): The DataFrame containing projected moment columns.
    Returns:
        pd.DataFrame: The DataFrame with the specified projected moment columns renamed.
    """
    return df.rename(
        columns={
            col: col.rsplit("_", 1)[0]
            for col in df.columns
            if col.endswith("_real") or col.endswith("_imag")
        }
    )
