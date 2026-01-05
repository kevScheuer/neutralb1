/**
 * @file fit_utils.h
 * @author Kevin Scheuer
 * @brief Utility functions for handling AmpTools fit results
 *
 */

#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "IUAmpTools/FitResults.h"
#include "TH1.h"

#include "neutralb1/AmplitudeParser.h"

/**
 * @brief Calculate the parity of the omega pi0 system
 *
 * @param[in] L Angular momenta between the omega and pi0
 * @return int The calculated parity
 */
int calculate_system_parity(int L);

/**
 * @brief Find all amplitudes matching the AmplitudeParser object passed.
 *
 * @details
 * At minimum, the AmplitudeParser object should have the e, J, P, m, and L values, and
 * so this function will return all amplitudes in the FitResults that match those
 * quantum numbers. If the reaction and/or sum fields of the AmplitudeParser object
 * are also filled, then the search will be restricted to those as well.
 *
 * @param[in] amp_to_find AmplitudeParser object specifying the quantum numbers (and
 *  optionally reaction and sum) to search for
 * @param[in] results FitResults object containing the amplitudes to search through
 * @param[in] skip_background Whether to skip isotropic background amplitudes. Defaults
 *  to true.
 * @return std::vector<std::string> Full amplitude strings matching amp_to_find
 */
std::vector<std::string> find_matching_amplitudes(
    const AmplitudeParser &amp_to_find,
    const FitResults &results,
    const bool skip_background = true);

/**
 * @brief Get the production coefficient for the wave's quantum numbers.
 *
 * @details
 * Since we aren't explicitly looping over the parity values, we'll need to infer them
 * from the J and L values. This function is thus hard-coded for omega-pi production
 * processes for now, and is reaction and sum independent. Any fit with multiple
 * reactions and non-constrained sums could be subject to undefined behavior.
 *
 * @param[in] e reflectivity of the first wave
 * @param[in] J total angular momentum of the first wave
 * @param[in] m m-projection of the first wave
 * @param[in] L orbital angular momentum of the first wave
 * @param[in] e_conj reflectivity of the conjugate wave
 * @param[in] J_conj total angular momentum of the conjugate wave
 * @param[in] m_conj m-projection of the conjugate wave
 * @param[in] L_conj orbital angular momentum of the conjugate wave
 * @param[in] reaction The reaction string (for polarization orientation)
 * @param[in] results The fit results containing the necessary data.
 * @param[in] acceptance_corrected Whether to use acceptance-corrected production coefficients. Defaults to true.
 * @return The production coefficient if found, or 0 if not found.
 */
complex<double> get_production_coefficient_pair(
    int e, int J, int m, int L,
    int e_conj, int J_conj, int m_conj, int L_conj,
    const std::string &reaction, const FitResults &results,
    bool acceptance_corrected = true);

/**
 * @brief Calculate the intensity of the fit results.
 *
 * Useful for comparing to H0_0000 moment, or the AmpTools reported intensity, to check
 * that \ref get_production_coefficients "the production coefficient calculator" is
 * working properly.
 *
 * @param[in] results The fit results to calculate the intensity for.
 * @param[in] acceptance_corrected Whether to return the acceptance-corrected intensity. Defaults to true.
 * @return complex<double> The calculated intensity.
 */
double calculate_intensity(const FitResults &results, bool acceptance_corrected = true);

/**
 * @brief Find the maximum spin J value from the fit results.
 * @param[in] results The fit results to search through.
 * @return int The maximum J value found in the fit results.
 */
int find_max_J(const FitResults &results);

/**
 * @brief Get the bin width of a 1D histogram
 *
 * @param[in] h Pointer to the histogram
 * @return double The bin width of the histogram
 */
inline double get_bin_width(TH1F *h)
{
    return (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) / h->GetNbinsX();
}

/**
 * @brief Join the keys of a color map into a single TString separated by a delimiter
 *
 * @param m map of cut TStrings to color integers
 * @param delimiter delimiter to separate keys, default is ","
 * @return TString joined keys
 */
TString join_keys(const std::map<TString, Int_t> &m, const TString &delimiter = ",");

/**
 * @brief Structure to hold vector-pseudoscalar moment quantum numbers
 *
 * @details The five quantum numbers needed to specify a moment:
 * - alpha: indexes the intensity term
 *  0: unpolarized, 1: polarized cos(2*Phi), 2: polarized sin(2*Phi)
 * - Jv: related to spin of the vector meson
 * - Lambda: related to helicity of the vector meson
 * - J: total angular momentum of the moment
 * - M: spin projection of the moment
 */
struct Moment
{
    int alpha;
    int Jv;
    int Lambda;
    int J;
    int M;

    std::string name() const
    {
        return "H" +
               std::to_string(alpha) + "_" +
               std::to_string(Jv) +
               std::to_string(Lambda) +
               std::to_string(J) +
               std::to_string(M);
    }

    bool operator<(const Moment &other) const
    {
        if (alpha != other.alpha)
            return alpha < other.alpha;
        if (Jv != other.Jv)
            return Jv < other.Jv;
        if (Lambda != other.Lambda)
            return Lambda < other.Lambda;
        if (J != other.J)
            return J < other.J;
        return M < other.M;
    }
    bool operator==(const Moment &other) const
    {
        return alpha == other.alpha && Jv == other.Jv &&
               Lambda == other.Lambda && J == other.J &&
               M == other.M;
    }
};

/**
 * @brief Initialize the set of moments based on the fit results.
 *
 * @param[in] results The fit results containing the necessary data.
 * @return std::vector<Moment> A vector of all possible moments constructed from the
 *  fit results.
 */
std::vector<Moment> initialize_moments(const FitResults &results);

#endif // FIT_UTILS_H