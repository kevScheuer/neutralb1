#include "file_utils.h"
#include <fstream>
#include <iostream>

std::vector<std::string> read_file_list(const std::string& file_path) {
    std::vector<std::string> file_vector;
    std::ifstream infile(file_path);
    
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return file_vector; // returns empty vector
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines
        if (!line.empty()) {
            file_vector.push_back(line);
        }
    }
    
    return file_vector;
}
