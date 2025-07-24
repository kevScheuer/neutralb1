"""move here
- bootstrap testing
- statistical analysis

make sure bootstrap methods are grouped by the file directory
"""

# TODO: complete this function. Uses old variables from old Plotter class
# def normality_test(
#     self,
#     columns: Optional[List[str]] = None,
#     method: Optional[Literal["None", "normal", "log_normal"]] = None,
# ) -> None:
#     """
#     Shapiro-wilkes for parameters and amplitudes. Log-transform the amps near zero.
#     Use other method for circular (phase) data.

#     """
#     if self.bootstrap_df is None:
#         raise ValueError("Bootstrap df was not defined on instantiation")

#     # PLAN: loop over the fit indices of the bootstrap df and calculate the test.
#     # Those that don't pass, or whose bootstrap mean is too far from the fit
#     # result's are plotted in big pdf. There should be indicators for the fit mean
#     # value on the QQ plot, and whether it failed due to p value or mean

#     fit_indices = self.bootstrap_df["fit_index"].unique()

#     # add the amplitudes and phase differences to the columns to be tested
#     columns = columns if columns is not None else []
#     columns.extend(self.coherent_sums["eJPmL"])
#     columns.extend(set(self.phase_differences.values()))

#     failed_dict = {}  # to store {fit_index: {column : [mean, p-value]}}
#     for fit_index in fit_indices:
#         for col in columns:
#             if col in set(self.phase_differences.values()):
#                 # do circular data test. Take absolute value because we have no way
#                 # to distinguish positive vs negative values
#                 pass

#             elif any(col in sublist for sublist in self.coherent_sums.values()):
#                 # do normality test on amplitudes
#                 # log-transform the amps near zero
#                 # type: ignore or type annotation for pylance
#                 amp_values = self.bootstrap_df.loc[
#                     self.bootstrap_df["fit_index"] == fit_index
#                 ][col]

#                 # log-transform if col is near zero
#                 # use rel distance to minimum test to determine if near zero
#                 # TODO: finish this

#                 # perform normality test
#                 stat, p_value = scipy.stats.shapiro(amp_values)
#                 if p_value < 0.05:
#                     failed_dict.setdefault(fit_index, {})[col] = [
#                         amp_values.mean(),
#                         p_value,
#                     ]

#             else:
#                 # do normality test as usual
#                 values = self.bootstrap_df.loc[
#                     self.bootstrap_df["fit_index"] == fit_index
#                 ][col]

#                 # perform normality test
#                 stat, p_value = scipy.stats.shapiro(values)
#                 if p_value < 0.05:
#                     failed_dict.setdefault(fit_index, {})[col] = [
#                         values.mean(),
#                         p_value,
#                     ]

#     # PLAN: here we'll loop over the failed ones and plot probdists

#     pass
