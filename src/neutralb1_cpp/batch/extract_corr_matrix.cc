/* Obtain the correlation matrix from AmpTools .fit results and save to a csv file

This script must be written in ROOT/c++ in order to interface with the AmpTools library
and access the FitResults class, which will be how we extract the correlation matrix.
The csv file produced will contain all correlation values, with the labelled columns and
rows corresponding to the parameters

This is essentially a copy of extract_cov_matrix.cc, but written in a separate file to
prevent the user from having to manually calculate the correlations by accessing all the
parameter errors.

NOTE: This script assumes that 4 coherent sums, 2 for each reflectivity, are present,
and that the reflectivities are constrained to each other. Thus we don't print 1/2 of
the reflectivity dependent parameters, as they're information is fixed to the other
half. If we left them in, we would have lots of redundant info and artificially highly
covariant parameters. If you have custom fit models that don't follow this rule, then
remove any lines that "continue" when the "Imag" substring is found.
*/

#include <cstring>
#include <fstream> // for writing csv
#include <iostream>
#include <sstream> // for std::istringstream
#include <string>
#include <vector>

#include "IUAmpTools/FitResults.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <file_list.txt> <output.csv>\n";
        return 1;
    }

    std::string file_path = argv[1];
    std::string csv_name = argv[2];

    // file path is a text file with a list of AmpTools output files, each on a newline
    std::vector<std::string> file_vector;
    std::ifstream infile(file_path);
    std::string line;
    while (std::getline(infile, line))
    {
        file_vector.push_back(line);
    }

    std::ofstream csv_file;
    csv_file.open(csv_name);

    // Collect all rows in a stringstream to minimize I/O operations
    std::stringstream csv_data;

    // ==== BEGIN FILE ITERATION ====
    // Iterate over each file, and label each correlation matrix block with the file
    // name
    bool is_header_written = false; // only write header once
    long unsigned int num_params = 0; // store number of parameters to ensure other files match
    for (const std::string &file : file_vector)
    {
        std::cout << "Analyzing File: " << file << "\n";
        FitResults results(file);
        if (!results.valid())
        {
            std::cout << "Invalid fit results in file: " << file << "\n";
            continue;
        }
        // write the header row
        if (!is_header_written)
        {
            is_header_written = true;
            csv_data << "file" << ",";
            csv_data << "parameter";
            
            num_params = results.parNameList().size();
            for (const auto &par : results.parNameList())
            {
                // skip parameters constrained across coherent sums. We only need the
                // "Real" coherent sums, as the "Imag" ones are constrained to them.
                // (ImagNegSign <-> RealPosSign) and (ImagPosSign <-> RealNegSign).
                if (par.find("Imag") != std::string::npos)
                    continue;
                csv_data << "," << par;
            }
            csv_data << "\n";
        } // end header row

        // get the covariance matrix
        const std::vector<std::vector<double>> cov_matrix = results.errorMatrix();

        // check that the covariance matrix matches the number of parameters
        if (cov_matrix.size() != results.parNameList().size())
        {
            std::cout << "Error: Covariance matrix size does not match number of parameters. Exiting! \n";
            exit(1);
        }
        // check that this file has the same number of parameters as the first file
        if (num_params != results.parNameList().size())
        {
            std::cout << "Error: Number of parameters in file " << file
                      << " does not match those of " << file_vector[0] 
                      << " Exiting! \n";
            exit(1);
        }

        for (long unsigned int row = 0; row < cov_matrix.size(); row++)
        {
            const std::string row_par = results.parNameList()[row];
            const double row_par_error = results.parError(row_par);
            if (row_par.find("Imag") != std::string::npos)
                continue;

            // repeat file name to label each block of the correlation matrix, and
            // then write the parameter name in the index column
            csv_data << file << ",";
            csv_data << row_par;

            for (long unsigned int col = 0; col < cov_matrix[row].size(); col++)
            {
                const std::string col_par = results.parNameList()[col];
                const double col_par_error = results.parError(col_par);
                if (col_par.find("Imag") != std::string::npos)
                    continue;

                // check to ensure that our indexing is right by comparing the
                // covariance value accessed by index to the covariance accessed by
                // parameter name
                if (cov_matrix[row][col] != results.covariance(row_par, col_par))
                {
                    std::cout << "Error: Mismatch in covariance values between parameters "
                         << row_par << " and " << col_par << ". Exiting! \n";
                    exit(1);
                }

                // check that denominator is not zero to avoid division by zero
                if (row_par_error == 0.0 || col_par_error == 0.0)
                {
                    csv_data << ",0.0"; // set correlation to 0 if error is zero
                    continue;
                }
                double correlation = cov_matrix[row][col] / (row_par_error * col_par_error);
                // check for NaN or Inf values
                if (std::isnan(correlation) || std::isinf(correlation))
                {
                    std::cerr << "Error: Correlation is NaN or Inf between parameters "
                              << row_par << " and " << col_par << ". Exiting!\n";
                    throw std::runtime_error("Invalid correlation value detected.");
                }
                csv_data << "," << correlation;
            }
            csv_data << "\n";
        }
        // cov matrix written, move on to next file
    }

    // Write all collected data to the CSV file at once
    csv_file << csv_data.str();
    csv_file.close();

    return 0;
}