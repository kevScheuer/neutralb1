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
        proj_moments_df: Optional[pd.DataFrame] = None,
        randomized_df: Optional[pd.DataFrame] = None,
        randomized_proj_moments_df: Optional[pd.DataFrame] = None,
        bootstrap_df: Optional[pd.DataFrame] = None,
        bootstrap_proj_moments_df: Optional[pd.DataFrame] = None,
        truth_df: Optional[pd.DataFrame] = None,
        truth_proj_moments_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Initialize the factory with common data and utilities."""
        self.fit_df = fit_df
        self.data_df = data_df
        self.proj_moments_df = proj_moments_df
        self.randomized_df = randomized_df
        self.randomized_proj_moments_df = randomized_proj_moments_df
        self.bootstrap_proj_moments_df = bootstrap_proj_moments_df
        self.bootstrap_df = bootstrap_df
        self.truth_df = truth_df
        self.truth_proj_moments_df = truth_proj_moments_df

        # Set the matplotlib style for consistent plotting
        WORKSPACE_DIR = utils.get_workspace_dir()
        plt.style.use(f"{WORKSPACE_DIR}/config/neutralb1.mplstyle")

    @property
    def intensity(self):
        return IntensityPlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            proj_moments_df=self.proj_moments_df,
            randomized_df=self.randomized_df,
            randomized_proj_moments_df=self.randomized_proj_moments_df,
            bootstrap_df=self.bootstrap_df,
            bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
            truth_df=self.truth_df,
            truth_proj_moments_df=self.truth_proj_moments_df,
        )

    @property
    def phase(self):
        return PhasePlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            proj_moments_df=self.proj_moments_df,
            randomized_df=self.randomized_df,
            randomized_proj_moments_df=self.randomized_proj_moments_df,
            bootstrap_df=self.bootstrap_df,
            bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
            truth_df=self.truth_df,
            truth_proj_moments_df=self.truth_proj_moments_df,
        )

    @property
    def diagnostic(self):
        return DiagnosticPlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            proj_moments_df=self.proj_moments_df,
            randomized_df=self.randomized_df,
            randomized_proj_moments_df=self.randomized_proj_moments_df,
            bootstrap_df=self.bootstrap_df,
            bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
            truth_df=self.truth_df,
            truth_proj_moments_df=self.truth_proj_moments_df,
        )

    @property
    def randomized(self):
        return RandomizedPlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            proj_moments_df=self.proj_moments_df,
            randomized_df=self.randomized_df,
            randomized_proj_moments_df=self.randomized_proj_moments_df,
            bootstrap_df=self.bootstrap_df,
            bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
            truth_df=self.truth_df,
            truth_proj_moments_df=self.truth_proj_moments_df,
        )

    @property
    def bootstrap(self):
        return BootstrapPlotter(
            fit_df=self.fit_df,
            data_df=self.data_df,
            proj_moments_df=self.proj_moments_df,
            randomized_df=self.randomized_df,
            randomized_proj_moments_df=self.randomized_proj_moments_df,
            bootstrap_df=self.bootstrap_df,
            bootstrap_proj_moments_df=self.bootstrap_proj_moments_df,
            truth_df=self.truth_df,
            truth_proj_moments_df=self.truth_proj_moments_df,
        )
