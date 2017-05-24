#ifndef __SIGMA_IPOPT_RUN_H__
#define __SIGMA_IPOPT_RUN_H__

//=============================================================================
// sigma_ipopt_run.h
//   : Initialize and run IPOPT with options.
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
//#include "parse_config.h"
//#include "sigma_config.h"
//#include "sigma_core.h"

// run IPOPT
int runIpopt(
             const std::string & working_directory,    // (in) working directory
             const std::string & input_qmatrix,        // (in) qmatrix filename
             const int number_processes);              // (in) number of threads

#endif
