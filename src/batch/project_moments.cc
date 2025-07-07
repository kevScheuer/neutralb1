/*Project vector-pseudoscalar PWA fit results to moments

This script uses the conversion from partial wave complex values to "project" PWA fit
results into unique moments. The file take AmpTools .fit files as input, which contain
the results of a PWA fit. The moments are computed and saved to an output csv file. 

NOTE: This script assumes that the amplitudes are written in the vec-ps eJPmL format.
For example, the positive reflectivity, JP=1+, m=0, S-wave amplitude would be written
in the cfg file as [reaction]::RealNegSign::p1p0S. If you have a different format, then
you'll have to account for it in parse_amplitude()

TODO: Fill in
Usage: project_moments  
*/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "file_utils.h"


int main(int argc, char* argv[]) 
{
    std::string input_file = argv[1];
    std::string csv_name = argv[2];

    // input file is a text file with a list of .fit results, each on a newline.
    // load this into a vector using the utility function
    std::vector<std::string> file_vector = read_file_list(input_file);
    
    if (file_vector.empty()) {
        std::cerr << "Error: Could not read file list from " << input_file << std::endl;
        return 1;
    }
}