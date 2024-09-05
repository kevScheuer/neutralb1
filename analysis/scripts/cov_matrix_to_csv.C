/* Obtain the covariance matrix from an AmpTools .fit result file and save to a csv file

This script must be written in ROOT/c++ in order to interface with the AmpTools library
and access the FitResults class, which will be how we extract the covariance matrix.
The csv file produced will contain all covariance values, with the columns and rows
corresponding to the parameters
*/

#include "IUAmpTools/FitResults.h"

/*
    Args:
        file (string): AmpTools .fit file to extract matrix from
        csv_name (string, optional): Name of the csv file to be created. Defaults to
            'covariance.csv'
*/
void cov_matrix_to_csv(std::string file, std::string csv_name = "covariance.csv")
{
    FitResults results(file); // load fit result from file
    if (!results.valid())
    {
        cout << "Invalid fit results in file: " << file << ". Exiting \n";
        exit(1);
    }

    ofstream csv_file;
    csv_file.open(csv_name);

    // write header line
    csv_file << "parameters"; // this will be our csv "index" column
    for(auto par : results.parNameList())
    {
        // Parameters across sums are constrained such that 
        // (ImagNegSign <-> RealPosSign) and (ImagPosSign <-> RealNegSign).
        // Params from constrained sums will appear to be highly correlated and are
        // unnecessary to include. So any "Imag" param is skipped over
        if(par.find("Imag") != std::string::npos)
            continue;
        csv_file << "," << par;    
    }
        
    csv_file << "\n";
    
    std::vector<std::vector<double>> cov_matrix = results.errorMatrix();
    for(int row=0; row<cov_matrix.size(); row++) 
    {        
        std::string row_par = results.parNameList()[row];
        if(row_par.find("Imag") != std::string::npos)
            continue;

        // write parameter name in index column
        csv_file << row_par;

        for(int col=0; col<cov_matrix[row].size(); col++) 
        {
            std::string col_par = results.parNameList()[col];            
            if(col_par.find("Imag") != std::string::npos)
                continue;

            // check to ensure that right values are being called cov_matrix by 
            // explicitly requesting the covariance between the row and column params
            if(cov_matrix[row][col] != results.covariance(row_par, col_par))
            {
                cout << "Error: Mismatch in covariance values between parameters " 
                    << row_par << " and " << col_par << ". Exiting! \n";
                exit(1);
            }

            csv_file << "," << cov_matrix[row][col];
        }
        csv_file << "\n";
    }

    return;
}