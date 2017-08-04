#ifndef __IPOPT_NLP_H__
#define __IPOPT_NLP_H__

//=============================================================================
// ipopt_nlp.h
//   : This is the header file for solving model using IPOPT
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "IpTNLP.hpp"
#include <set>

using namespace Ipopt;

class IPOPT_NLP : public TNLP
{
public:
    // Default constructor
    IPOPT_NLP(const std::string & input_filepath, 
                     const int & number_threads);

    // Default destructor
    virtual ~IPOPT_NLP();

    // Method to return some info about the nlp
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style);

    // Method to return the bounds for my problem
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    // Method to return the starting point for the algorithm 
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda);

    // Method to return the objective value
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

    // Method to return the gradient of the objective
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    // Method to return the constraint residuals
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    // Method to return:
    //   1) The structure of the jacobian (if "values" is NULL)
    //   2) The values of the jacobian (if "values" is not NULL)
    //
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values);

    // Method to return:
    //   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    //   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    // 
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);

    // This method is called when the algorithm is complete so the TNLP can store/write the solution 
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);

private:
    //  IPOPT_NLP();
    IPOPT_NLP(const IPOPT_NLP&);
    IPOPT_NLP& operator=(const IPOPT_NLP&);
    //@}

    //unsigned int unique_number_aligned_reads;
    unsigned int unique_number_aligned_reads;

    // unique number of aligned genomes
    Index unique_number_aligned_genomes;

    // unique number of aligned genomes
    std::string output_filename;

    // vector of vector for aligned genome index
    std::vector< std::vector<unsigned int> > *aligned_genome_index_data;

    // vector of vector for aligned q-value 
    std::vector< std::vector<double> > *aligned_qvalue_data;

    // vector for genome index info
    std::vector<std::string> *genome_index_info;

    // for total number of reads
    std::vector<std::string> *reads_info;
    std::vector<double> *reads_weights;

    // for minimum_percentage_chance_genome
    double min_relative_abundance_genome;
};

#endif
