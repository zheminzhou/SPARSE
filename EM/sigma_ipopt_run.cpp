//=============================================================================
// sigma_ipopt_run.cpp
//   : Initialize IPOPT with options.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <omp.h>

#include "IpIpoptApplication.hpp"
#include "sigma_ipopt_nlp.h"
#include "utils.h"

//=============================================================================
// run IPOPT
//=============================================================================
int runIpopt(
             const std::string & working_directory, 
             const std::string & input_qmatrix, 
             const int number_threads)
{
    // [Step 1] prepare
    std::cout << "  ** Check the qmatrix model: Running -> " << std::flush;

    // initialize
    std::string qmatrix_filename, qmatrix_filebase, qmatrix_filepath;

    // get qmatrix filename and path
    if (input_qmatrix != "") 
        qmatrix_filename = input_qmatrix; // Utils::getFilename(input_qmatrix);
    else
        qmatrix_filename = "sigma_out.qmatrix.txt";
    qmatrix_filebase = Utils::getFilebase2(qmatrix_filename);
    qmatrix_filepath = qmatrix_filename;

    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    SmartPtr<TNLP> mynlp = new SIGMA_IPOPT_NLP(qmatrix_filepath, 
                                               qmatrix_filebase, 
                                               number_threads);

    // for ipopt output
    std::string ipopt_output = qmatrix_filebase + ".ipopt.txt";

    // [Step 1] done
    std::cout << "Done!" << std::endl;

    // [Step 2] align reads to each genome
    std::cout << "  ** Solve the model using IPOPT: Running -> " << std::flush;

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    SmartPtr<IpoptApplication> app = new IpoptApplication();

    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) {
        std::cout << "  ** The problem solved!" << std::endl << std::endl;
    }
    else {
        std::cout << "  ** The problem FAILED!" << std::endl << std::endl;
    }

    return (int) status;
}
