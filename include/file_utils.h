#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <vector>
#include <string>

/**
 * Reads a text file containing a list of file paths (one per line) and returns them as a vector
 * @param file_path Path to the text file containing the list of files
 * @return Vector of file paths as strings
 */
std::vector<std::string> read_file_list(const std::string &file_path);

#endif // FILE_UTILS_H
