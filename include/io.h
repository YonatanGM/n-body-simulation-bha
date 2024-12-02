#ifndef IO_H
#define IO_H

#include <vector>
#include <string>
#include "body.h"

// function to read bodies from a CSV file
bool readCSV(const std::string& filename, std::vector<Body>& bodies);

// function to save the state of bodies to a CSV file
void saveState(const std::string& vs_dir, int vs_counter, const std::vector<Body>& bodies);

// Function to convert orbital elements CSV to state vector CSV
bool convertOrbitalElementsToCSV(const std::string& inputFilename, const std::string& outputFilename);

// Function to combine two CSV files into one, removing duplicates based on the "name" column
void combineCSVFiles(const std::string& inputFile1, const std::string& inputFile2, const std::string& outputFile);

#endif // IO_H
