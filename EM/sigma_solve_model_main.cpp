//=============================================================================
// sigma_solve_model_main.cpp
//   : main cpp code to call IpOpt for solving non-linear problem
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <omp.h>

#include "IpIpoptApplication.hpp"
#include "sigma_ipopt_run.h"
#include "sigma_ipopt_nlp.h"
#include "utils.h"


using namespace Ipopt;


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name);


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name)
{
    std::cout << std::endl
              << "  [Usage]" << std::endl 
              << "    " << program_name << " [options] -w <working directory>" << std::endl
              << std::endl
              << "  [Inputs]" << std::endl
              << "    1. working directory (default: current running directory)" << std::endl
              << "      - if working_directory is not specified, the program will work in the current directory" << std::endl
              << "      - results will be generated in the working directory" << std::endl
              << std::endl
              << "  [Options]" << std::endl
              << "    -h/--help" << std::endl
              << "    -v/--version" << std::endl
              << "    -i/--input-qmatrix <string> # provide q-matrix filename directry" << std::endl
              << "    -t/--multi-threads <int>    # number of threads (default: 1)" << std::endl
              << std::endl
              << "  [Outputs]" << std::endl
              << "    sigma_out.gvector.txt" << std::endl
              << "    sigma_out.gvector.html" << std::endl
              << "    sigma_out.ipopt.txt" << std::endl
              << std::endl;
}


//=============================================================================
// Initialize arguments
//=============================================================================
void initializeArguments(
        int argc, char **argv,              // (in) argv
        std::string & working_directory,    // (out) working directory
        std::string & input_qmatrix,        // (out) input qmatrix
        int & number_threads)               // (out) number of multi-threads
{
    // initialize variables
    int i, number_threads_default;
    std::string  program_name, program_path, input_qmatrix_default;
    std::string working_directory_default;

    // grab command line arguments
    std::vector<std::string> arguments_vector;

    // push to arguments_vector
    while(argc--) {
        arguments_vector.push_back(*argv++);
    }

    // default option values
    working_directory_default = std::string(".") + Utils::getPathSeparator();
    input_qmatrix_default     = "sigma_out.qmatrix.txt";
    number_threads_default  = 1;

    // option values
    working_directory = working_directory_default;
    number_threads  = number_threads_default;
    input_qmatrix     = input_qmatrix_default;

    // get program name
    program_path = arguments_vector[0];
    program_name = Utils::getProgramName(program_path);

    // get argements
    for(i = 1; i <= (int)arguments_vector.size()-1; i++)
    {   
        // working directory
        if(arguments_vector[i] == "-w") {
            working_directory = arguments_vector[++i];
            if (*working_directory.rbegin() != Utils::getPathSeparator())
                working_directory += Utils::getPathSeparator();
        }
        // config path
        else if (arguments_vector[i] == "-i" || 
                 arguments_vector[i] == "--input-qmatrix") {
            input_qmatrix = arguments_vector[++i];
        }
        // multi-threads
        else if (arguments_vector[i] == "-t" || 
                 arguments_vector[i] == "--multi-threads") {
            std::stringstream(arguments_vector[++i]) >> number_threads;
            if (number_threads < 1) 
                Utils::exitWithError("*** Error: check -t option value.");
        }
        // help usage
        else if (arguments_vector[i] == "-h" ||
                 arguments_vector[i] == "--help") {
            usage(program_name);
            exit(0);
        }
        // unknown option
        else {
            std::cerr << "*** Error: Unknown option " << arguments_vector[i]
                << std::endl << std::endl;
            usage(program_name);
            exit(1);
        }
    }
}


//=============================================================================
// Main
//=============================================================================
int main(int argc, char** argv)
{

    // get program name
    std::string program_path = argv[0];
    std::string program_name = "SPARSED EM";

    // initialize variables
    int number_threads;
    double start_time, finish_time, elapsed_time;
    std::string working_directory, input_qmatrix;

    // initialize arguments
    initializeArguments(argc, argv, // (in)
        working_directory, input_qmatrix, number_threads); // (out)

    // record start time
    start_time = omp_get_wtime();

    // display work start and time record
    std::cout << std::endl
        << "********************************************************************************" << std::endl
        << Utils::currentDateTime() << " Beginning. " <<  std::endl;

    //=============================================================================
    // run IPOPT
    //=============================================================================
    runIpopt(working_directory, input_qmatrix, number_threads);

    // finish time
    finish_time = omp_get_wtime();
    elapsed_time = double(finish_time - start_time);

    // display elapsed time
    std::cout << Utils::currentDateTime() << " Ending " << program_name << std::endl
        << "Total Elapsed Time =  " << elapsed_time << " [seconds]" << std::endl
        << "********************************************************************************" 
        << std::endl << std::endl;
}
