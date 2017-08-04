#ifndef __IPOPT_RUN_H__
#define __IPOPT_RUN_H__

//=============================================================================
// ipopt_run.h
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
#include "ipopt_nlp.h"

// run IPOPT
int runIpopt(
             const std::string & input_qmatrix,        // (in) qmatrix filename
             const int number_processes);              // (in) number of threads

#endif
