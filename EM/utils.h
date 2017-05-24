#ifndef UTILS_H
#define UTILS_H

//=============================================================================
// utils.h
//   : header file for utils.cpp that includes useful common uility functions.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <functional> 
#include <iostream>
#include <locale>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>


//=============================================================================
// utils: utility functions for file (directory) path, filename,
//        base filename without extension, and format changing.
//=============================================================================
namespace Utils
{
    // Get path separator
    char getPathSeparator();

    // Check file exist
    bool isFileExist(const std::string & filename);

    // check executable file
    bool canExec(const std::string & filepath);

    // If file exists, then remove the file
    void ifFileExistRemove(const std::string & filepath);

    // Get program name only from program path
    std::string getProgramName(const std::string & program_path);

    // Get program directory from program path
    std::string getProgramDir(const std::string & program_path);

    // If error, then print error and exit program 
    void exitWithError(const std::string & error);

    // Method to convert string to int 
    int stringToInt(const std::string & input_string);

    // Method to convert string to unsigned int
    unsigned int stringToUnsignedInt(const std::string & input_string);

    // Method to convert string to double
    double stringToDouble(const std::string & input_string);

    // Method to convert int to string 
    std::string intToString(const int & input_value);

    // Method to convert unsigned int to string 
    std::string unsignedIntToString(const unsigned int & input_value);

    // Method to convert double to string 
    std::string doubleToString(const double & input_value);

    // Check directory 
    bool isDirectory(const std::string & directory);

    // Check fasta format 
    bool isFastaFormat(const std::string & input_string);

    // Check fastq format 
    bool isFastqFormat(const std::string & input_string);

    // Get file name only from file path 
    std::string getFilename(const std::string & filepath);

    // Get filebase from filename without extension 
    std::string getFilebase(const std::string & filename);

    // Get filebase from filename without extension 
    std::string getFilebase2(const std::string & filename);

    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss 
    std::string currentDateTime();

    // String vector to comma separated string 
    std::string vectorToCommaString(const std::vector<std::string> & input_vector);

    // Comma separated string to vector
    std::vector<std::string> commaStringToVector(const std::string & input_string);

    /*
    // trim from start
    std::string ltrim(std::string &s);

    // trim from end
    std::string rtrim(std::string &s);

    // trim from both ends
    std::string trim(std::string &s);
    */

};

#endif
