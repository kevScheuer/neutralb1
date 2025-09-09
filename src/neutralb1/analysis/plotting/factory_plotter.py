from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

import neutralb1.utils as utils
from neutralb1.analysis.plotting.bootstrap_plotter import BootstrapPlotter
from neutralb1.analysis.plotting.diagnostic_plotter import DiagnosticPlotter
from neutralb1.analysis.plotting.intensity_plotter import IntensityPlotter
from neutralb1.analysis.plotting.phase_plotter import PhasePlotter
from neutralb1.analysis.plotting.randomized_plotter import RandomizedPlotter


class FactoryPlotter:
    """Factory class that interfaces with all sub-plotters."""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        random_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Initialize the factory with common data and utilities."""
        self.fit_df = fit_df
        self.data_df = data_df
        self.random_df = random_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df

        # Set the matplotlib style for consistent plotting
        WORKSPACE_DIR = utils.get_workspace_dir()
        plt.style.use(f"{WORKSPACE_DIR}/config/neutralb1.mplstyle")

    def intensity(self):
        return IntensityPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )

    def phase(self):
        return PhasePlotter(self.fit_df, self.data_df, self.bootstrap_df, self.truth_df)

    def diagnostic(self):
        return DiagnosticPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )

    def random(self):
        return RandomizedPlotter(self.fit_df, self.data_df, self.random_df)

    def bootstrap(self):
        return BootstrapPlotter(
            self.fit_df, self.data_df, self.bootstrap_df, self.truth_df
        )
