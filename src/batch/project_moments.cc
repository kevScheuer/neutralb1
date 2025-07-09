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



// forward declarations
std::map<std::string, double> calculate_moments(const FitResults &results);

int main(int argc, char *argv[])
{
    std::string input_file = argv[1];
    std::string csv_name = argv[2];

    // TODO: Add a help message and argc check

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
    std::map<std::string, double> moments;

    // open csv file for writing
    std::ofstream csv_file;
    csv_file.open(csv_name);

    // Collect all rows in a stringstream to minimize I/O operations
    std::stringstream csv_data;
    bool is_header_written = false;

    // ==== BEGIN FILE ITERATION ====
    // Iterate over each file, and add their results as a row in the csv
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
        moments.clear();

        // calculate the moments from the fit results
        moments = calculate_moments(results);

    } // end of file iteration
}

std::map<std::string, double> calculate_moments(const FitResults &results)
{
    std::map<std::string, double> moments;

    return moments;
}

int find_max_m(const FitResults &results)
{
    int max_m = 0;
    for (const auto &reaction : results.reactionList())
    {
        for (const std::string &amplitude : results.ampList(reaction))
        {
            int m_value = std::stoi(parse_amplitude(amplitude).m);
            if (m_value > max_m)
            {
                max_m = m_value;
            }
        }
    }

    return max_m;
}