//=============================================================================
// utils.cpp
//   : utils.cpp includes useful common uility functions.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"


//=============================================================================
// utils: utility functions for file (directory) path, filename,
//        base filename without extension, and format changing.
//=============================================================================
namespace Utils
{

    // Get path separator 
    char getPathSeparator()
    {
    #if _WIN32
        return '\\' ;
    #else
        return '/' ;
    #endif
    }

    // Check file exist
    bool isFileExist(const std::string & filename)
    {
        // infile
        std::ifstream infile(filename.c_str());

        // return 
        return infile.good();
    }


    // Check executable file
    bool canExec(const std::string & filepath)
    {
        if (!access(filepath.c_str(), X_OK))
            return true;
        else
            return false;
        return false;
    }


    // If file exists, then remove the file
    void ifFileExistRemove(const std::string & filepath)
    {
        // if file exists,
        if (isFileExist(filepath))
        {
            // then remove if tere is no error
            if (remove(filepath.c_str()) != 0)
                exitWithError("*** Error: Failed to remove: " + filepath);
        }
    }


    // Get program name only from program path
    std::string getProgramName(const std::string & program_path)
    {
        std::string program_name;

        // position
        size_t program_position = program_path.find_last_of(getPathSeparator());

        // find path separator
        if (program_position != std::string::npos)
            program_name.assign(program_path.begin()+program_position+1,program_path.end());
        // couldn't find path separator
        else
            program_name = program_path;

        // return 
        return program_name;
    }

    // Get program directory from program path
    std::string getProgramDir(const std::string & program_path)
    {
        std::string program_dir;

        // position
        size_t program_position = program_path.find_last_of(getPathSeparator());

        // find path separator
        if (program_position != std::string::npos)
            program_dir.assign(program_path.begin(),program_path.begin()+program_position+1);
        // couldn't find path separator
        else
            program_dir = "." + getPathSeparator();

        // return 
        return program_dir;
    }

    // Exit program when error happens 
    void exitWithError(const std::string &error) 
    {
        std::cerr << std::endl << "  " << error << std::endl;
        exit(EXIT_FAILURE);
    }

    // Method to convert string to unsigned int
    int stringToInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        int x;
        if (!(input_string_stream >> x)) {
            exitWithError("*** Error: string was not converted to unsigned int correctly!");
        }   

        // return 
        return x;
    } 

    // Method to convert string to unsigned int
    unsigned int stringToUnsignedInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        unsigned int x;
        if (!(input_string_stream >> x)) {
            exitWithError("*** Error: string was not converted to unsigned int correctly!");
        }   

        // return 
        return x;
    } 

    // Method to convert string to double 
    double stringToDouble( const std::string& input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        double x;
        if (!(input_string_stream >> x)) {
            exitWithError("*** Error: string was not converted to double correctly!");
        }   

        // return 
        return x;
    } 

    // Method to convert int to string 
    std::string intToString( const int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    // Method to convert unsigned int to string
    std::string unsignedIntToString( const unsigned int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    // Method to convert double to string 
    std::string doubleToString( const double & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    // Check directory 
    bool isDirectory( const std::string & directory )
    {
        struct stat sb;
        if (stat(directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode) )
            return true;
        else
            return false;
    }

    // Check fasta format 
    bool isFastaFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fasta");
        file_format.push_back("fa");
        file_format.push_back("fas");
        file_format.push_back("fna");
        file_format.push_back("ffn");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);

        // transform to lower characters
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                return true;
        }
        // return 
        return false;
    }

    // Check fastq format 
    bool isFastqFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fastq");
        file_format.push_back("fq");
        file_format.push_back("faq");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                // return 
                return true;
        }
        // return 
        return false;
    }

    // Get file name only from file path 
    std::string getFilename(const std::string & filepath)
    {
        std::string filename;

        size_t pos = filepath.find_last_of(getPathSeparator());
        if (pos != std::string::npos)
            filename.assign(filepath.begin()+pos+1,filepath.end());
        else
            filename = filepath;

        // return 
        return filename;
    }

    // Get filebase from filename without extension 
    std::string getFilebase(const std::string & filename)
    {
        std::string filebase;

        size_t pos = filename.find_last_of(".");
        if (pos != std::string::npos)
            filebase.assign(filename.begin(),filename.begin()+pos);
        else
            filebase = filename;

        // return 
        return filebase;
    }

    // Get filebase from filename i.e., AAA.CC.X -> AAA 
    std::string getFilebase2(const std::string & filename)
    {
        std::string filebase;

        // parse base filename
        std::istringstream filename_stream(filename);
        std::vector<std::string> fields_vector;
        for (std::string field; getline(filename_stream, field, '.'); ) { 
           fields_vector.push_back(field);
        }   
        if (fields_vector.size() > 2) {
            filebase = ""; 
            for(unsigned int i = 0; i < fields_vector.size() - 2; i++) {
                filebase += fields_vector[i] + '.';
            }   
            filebase = filebase.substr(0, filebase.size()-1);
        }   
        else {
            filebase = filename;
        }

        // return 
        return filebase;
    }

    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss 
    std::string currentDateTime() 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "[%Y-%m-%d.%X]", &tstruct);

        // return 
        return buf;
    }

    // string vector to comma separated string 
    std::string vectorToCommaString(const std::vector<std::string> & input_vector)
    {
        // variables
        std::stringstream ss;
        std::string comma_string;

        // loop list
        for (unsigned int i=0; i<input_vector.size(); i++ ) {
            if (i != 0)
                ss << ",";
            ss << input_vector[i];
        }
  
        // comma_string
        comma_string = ss.str();

        // return 
        return comma_string;
    }

    // comma separated string to vector
    std::vector<std::string> commaStringToVector(const std::string & input_string)
    {
        // variables
        std::vector<std::string> elems;
        std::stringstream ss(input_string);
        std::string item;

        // loop 
        while (std::getline(ss, item, ',')) 
        {
            elems.push_back(item);
        }
  
        // return 
        return elems;
    }

    /*
    // trim from start
    std::string ltrim(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            return s;
    }

    // trim from end
    std::string rtrim(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
            return s;
    }

    // trim from both ends
    std::string trim(std::string &s) {
            return ltrim(rtrim(s));
    }
    */
}

