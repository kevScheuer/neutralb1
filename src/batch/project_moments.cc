/*Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to "project" PWA fit
results into unique moments. The file take AmpTools .fit files as input, which contain
the results of a PWA fit. The moments are computed and saved to an output csv file.

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S. If you have a different format, then
you'll have to account for it in parse_amplitude()

TODO: Fill in
TODO: Print a warning to the user if the fit results do not contain the same set of
    amplitudes. Have default behavior fill in the missing amplitudes with 0.
TODO: The current way is to get the production coefficient first, store them, then later
    fill them in when calculating the sdme. It might be easier to simply store the
    amplitude name then call the scaled production coefficient later?
Usage: project_moments
*/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <complex>
#include <map>
#include <unordered_set>

#include "IUAmpTools/FitResults.h"
#include "file_utils.h"
#include "amp_utils.h"

struct Moment
{
    int alpha;
    int Jv;
    int Lambda;
    int J;
    int M;

    std::string name() const
    {
        return "H" + std::to_string(alpha) + "_" + std::to_string(Jv) + std::to_string(Lambda) + std::to_string(J) + std::to_string(M);
    }
};

// forward declarations
std::vector<Moment> initialize_moments(const FitResults &results);
complex<double> calculate_moment(const Moment &moment, const FitResults &results);
int find_max_m(const FitResults &results);
int find_max_J(const FitResults &results);

int main(int argc, char *argv[])
{
    // Check if we have the required arguments
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_csv>" << std::endl;
        std::cerr << "  input_file: Text file containing list of .fit files" << std::endl;
        std::cerr << "  output_csv: Name of the output CSV file" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string csv_name = argv[2];

    // input file is a text file with a list of .fit results, each on a newline.
    // load this into a vector using the utility function
    std::vector<std::string> file_vector = read_file_list(input_file);

    if (file_vector.empty())
    {
        std::cerr << "Error: Could not read file list from " << input_file << std::endl;
        return 1;
    }

    // TODO: The following code should read a .fit from the vector, extract the moments,
    // and write them to the csv file.

    // initialize the map of moment names to their values
    std::map<std::string, complex<double>> moment_results;
    std::vector<Moment> moments; // vector of all moments to be calculated

    // open csv file for writing (force overwrite if file exists)
    std::ofstream csv_file(csv_name, std::ios::out | std::ios::trunc);

    // Collect all rows in a stringstream to minimize I/O operations
    std::stringstream csv_data;

    // ==== BEGIN FILE ITERATION ====
    for (const std::string &file : file_vector)
    {
        std::cout << "Analyzing File: " << file << "\n";
        FitResults results(file);
        if (!results.valid())
        {
            std::cout << "Invalid fit results in file: " << file << "\n";
            continue;
        }

        // before getting this file's info, clear the results from the last file
        moment_results.clear();
        moments.clear();

        // initialize the set of moments we can have from this file's waveset
        moments = initialize_moments(results);

        for (const Moment &moment : moments)
        {
            // calculate the value for this moment and save it to the map
            moment_results[moment.name()] = calculate_moment(moment, results);
        }

    } // end of file iteration
}

/**
 * @brief Initialize the set of moments based on the fit results.
 *
 * @param results The fit results containing the necessary data.
 * @return std::vector<Moment> A vector of all possible moments constructed from the
 *  fit results.
 */
std::vector<Moment> initialize_moments(const FitResults &results)
{
    std::vector<Moment> moments;
    int max_J = find_max_J(results);
    int max_m = find_max_m(results);

    // prepare moment quantum numbers for the moments
    std::vector<int> alpha_vector = {0, 1, 2};
    std::vector<int> Jv_vector = {0, 2}; // CGs coefficient ensure Jv=1 is always 0
    // moments with negative Lambda values are proportional to positive ones
    std::vector<int> Lambda_vector = {0, 1, 2};
    std::vector<int> J_vector;
    for (int J = 0; J <= max_J + 1; ++J) // J defined to be >= 0
    {
        J_vector.push_back(J);
    }
    std::vector<int> M_vector; // similar to lambda, (-) M values ar proportional to (+)
    for (int m = 0; m <= max_m + 1; ++m)
    {
        M_vector.push_back(m);
    }
    for (int alpha : alpha_vector)
    {
        for (int Jv : Jv_vector)
        {
            for (int Lambda : Lambda_vector)
            {
                for (int J : J_vector)
                {
                    for (int M : M_vector)
                    {
                        // FILTER non-physical moments
                        // Wigner D functions for these are always 0
                        if (M > J || Lambda > J || Lambda > Jv)
                            continue;
                        // Wishart seciton 5.10.2 proves H2(Jv,0,J,0) = 0 for any Jv, J
                        if (alpha == 2 && Lambda == 0 && M == 0)
                            continue;

                        // create a moment from these quantum numbers
                        Moment moment;
                        moment.alpha = alpha;
                        moment.Jv = Jv;
                        moment.Lambda = Lambda;
                        moment.J = J;
                        moment.M = M;

                        moments.push_back(moment);
                    }
                }
            }
        }
    }
    return moments;
}

/**
 * @brief Calculate the value of a moment based on its quantum numbers and fit results.
 *
 * @param moment The moment for which to calculate the value.
 * @param results The fit results containing the necessary data.
 * @return complex<double> The calculated value of the moment.
 */
complex<double> calculate_moment(const Moment &moment, const FitResults &results)
{
    complex<double> moment_value = 0.0;

    // for-loops below are done to best match the mathematical definition
    //  since moment.J = max(waveset J) + 1, we just need to iterate from -moment.J+1 to moment.J - 1
    for (int Ji = 0; Ji < moment.J; ++Ji)
    {
        for (int li = 0; li <= Ji; ++li)
        {
            for (int Jj = 0; Jj < moment.J; ++Jj)
            {
                for (int lj = 0; lj <= Jj; ++lj)
                {
                    for (int mi = -Ji; mi <= Ji; ++mi)
                    {
                        for (int mj = -Jj; mj <= Jj; ++mj)
                        {
                        }
                    }
                }
            }
        }
    }
}

// TODO: This should take in each group of Jlm pairs and calculate the SDME
//  I dont' like how the previous process_waves hid the conjugation and signs within it,
//  so I might do it all just within this function
complex<double> calculate_SDME(int alpha)
{
}

int find_max_m(const FitResults &results)
{
    int max_m = 0;
    for (const auto &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            if (amplitude.find("isotropic") != std::string::npos ||
                amplitude.find("Background") != std::string::npos)
                continue; // skip isotropic and background amplitudes

            int m_value = parse_amplitude(amplitude).get_m_int();
            if (m_value > max_m)
            {
                max_m = m_value;
            }
        }
    }

    return max_m;
}

int find_max_J(const FitResults &results)
{
    int max_J = 0;
    for (const auto &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            if (amplitude.find("isotropic") != std::string::npos ||
                amplitude.find("Background") != std::string::npos)
                continue; // skip isotropic and background amplitudes

            int J_value = parse_amplitude(amplitude).get_J_int();
            if (J_value > max_J)
            {
                max_J = J_value;
            }
        }
    }

    return max_J;
}