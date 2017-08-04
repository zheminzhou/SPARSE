//=============================================================================
// solve_model_main.cpp
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
#include "ipopt_run.h"
#include "ipopt_nlp.h"

using namespace Ipopt;



//=============================================================================
// Initialize arguments
//=============================================================================
void initializeArguments(
        int argc, char **argv,              // (in) argv
        std::string & input_qmatrix,        // (out) input qmatrix
        int & number_threads)               // (out) number of multi-threads
{
    // initialize variables
    int i, number_threads_default;
    std::string input_qmatrix_default;

    // grab command line arguments
    std::vector<std::string> arguments_vector;

    // push to arguments_vector
    while(argc--) {
        arguments_vector.push_back(*argv++);
    }

    // default option values
    input_qmatrix_default     = "sparse.qmatrix";
    number_threads_default  = 1;

    // option values
    number_threads  = number_threads_default;
    input_qmatrix     = input_qmatrix_default;

    // get argements
    for(i = 1; i <= (int)arguments_vector.size()-1; i++)
    {   
        // config path
        if (arguments_vector[i] == "-i" || 
                 arguments_vector[i] == "--input-qmatrix") {
            input_qmatrix = arguments_vector[++i];
        }
        // multi-threads
        else if (arguments_vector[i] == "-t" || 
                 arguments_vector[i] == "--multi-threads") {
            number_threads = std::stoul(arguments_vector[++i], nullptr);
            if (number_threads < 1) {
				std::cerr << std::endl << "*** Error: check -t option value." << std::endl;
				exit(EXIT_FAILURE);
            }
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
    std::string program_name = "sparse-solve-model";

    // initialize variables
    int number_threads;
    double start_time, finish_time, elapsed_time;
    std::string input_qmatrix;

    // initialize arguments
    initializeArguments(argc, argv, // (in)
        input_qmatrix, number_threads); // (out)

    // record start time
    start_time = omp_get_wtime();

    // display work start and time record
    std::cout << std::endl
        << "********************************************************************************" << std::endl
        << " Beginning. " <<  std::endl;

    //=============================================================================
    // run IPOPT
    //=============================================================================
    runIpopt(input_qmatrix, number_threads);

    // finish time
    finish_time = omp_get_wtime();
    elapsed_time = double(finish_time - start_time);

    // display elapsed time
    std::cout << " Ending " << program_name << std::endl
        << "Total Elapsed Time =  " << elapsed_time << " [seconds]" << std::endl
        << "********************************************************************************" 
        << std::endl << std::endl;
}
