from typing import Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import pandas as pd

from neutralb1.analysis.plotting.base_plotter import BasePWAPlotter


class RandomizedPlotter(BasePWAPlotter):
    """Handles plots that rely on randomized parameter fitting"""

    def __init__(
        self,
        fit_df: pd.DataFrame,
        data_df: pd.DataFrame,
        randomized_df: Optional[pd.DataFrame] = None,
        randomized_proj_moments_df: Optional[pd.DataFrame] = None,
    ) -> None:
        """Ensure random DataFrame is not None and initialize the plotter."""
        super().__init__(
            fit_df=fit_df,
            data_df=data_df,
            randomized_df=randomized_df,
            randomized_proj_moments_df=randomized_proj_moments_df,
        )
        if self.randomized_df is None:
            raise ValueError("RandomPlotter requires a non-None randomized_df.")

    def chi2_likelihood(self) -> matplotlib.axes.Axes:
        fig, ax = plt.subplots()

        # calculate chi2 / ndf based off intensities and phase differences between
        #   best fit and random samples
        # plot against Delta(-2ln(L)). Idea is that best fit is in bottom left corner,
        #   and other rand fits with different chi2 should also have different
        #   likelihoods
        # color code points corresponding to eMatrixStatus

        return ax
