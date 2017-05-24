//=============================================================================
// sigma_ipopt_nlp.cpp
//   : cpp code for solving model using IPOPT
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "sigma_ipopt_nlp.h"
#include "utils.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <ctime>
#include "omp.h"

using namespace Ipopt;


//=============================================================================
// Constructor 
//=============================================================================
SIGMA_IPOPT_NLP::SIGMA_IPOPT_NLP(const std::string & input_filepath, 
                                 const std::string & base_filename,
                                 const int & number_threads)
{
    // set number of threads
    omp_set_num_threads(number_threads);

    // define delimited characters
    char const tab_delimited = '\t';
    char const equal_delimited = '=';

    // typedef vector of vector
    typedef std::vector< std::vector<unsigned int> > vector_of_vector_unsigned_int;
    typedef std::vector< std::vector<double> > vector_of_vector_double;

    // allocate
    aligned_genome_index_data = new vector_of_vector_unsigned_int;
    aligned_qvalue_data = new vector_of_vector_double;
    genome_index_info = new std::vector<std::string>;
	reads_info = new std::vector<std::string>;
	reads_weights = new std::vector<double>;

    // open input file -> input_file
    std::ifstream input_file (input_filepath.c_str());

    // output filename
    output_filename = base_filename + ".gvector.txt";

    // max genome index
    unsigned int max_genome_index = 0;

    // get q-matrix data, and save to aligned_genome_index_data and aligned_qvalue_data
    for (std::string line; getline (input_file, line); ) 
    {
        // save genome index info to vector
		if (line[0] == '+') 
        {
            (*reads_info).push_back(line);
        }
        else if (line[0] == '@') 
        {
            (*genome_index_info).push_back(line);
        }
        else if (line[0] == '*') 
        {
            // line stringstream
            std::istringstream line_stream(line);

            // value_type for data
            (*aligned_genome_index_data).push_back(vector_of_vector_unsigned_int::value_type());
            (*aligned_qvalue_data).push_back(vector_of_vector_double::value_type());

            // loop token by tab-delimited and save to vector
            std::vector<std::string> fields_vector;
            for (std::string field; getline(line_stream, field, tab_delimited); ) {
                fields_vector.push_back(field);
            }

			(*reads_weights).push_back(Utils::stringToDouble(fields_vector[2]));
            // do not load the first column (ReadID) that will be used later
            for(unsigned int i = 3; i < fields_vector.size(); i++) {
                // field has pair "genome_index=q_value"
                std::string field=fields_vector[i];
                std::istringstream field_stream(field);

                // get genome_index and push to vector of vector
                unsigned int genome_index;
                std::string genome_index_string;
                getline(field_stream, genome_index_string, equal_delimited);
                genome_index = Utils::stringToUnsignedInt(genome_index_string);
                if (genome_index > max_genome_index) {
                    max_genome_index = genome_index;
                }
                (*aligned_genome_index_data).back().push_back(genome_index);
           
                // get genome_index and push to vector of vector
                double q_value;
                std::string q_value_string;
                getline(field_stream, q_value_string);
                q_value = Utils::stringToDouble(q_value_string);
                (*aligned_qvalue_data).back().push_back(q_value);
            }
            fields_vector.clear();
        }
    }

    // close input_file
    input_file.close();

    // get unique_number_aligned_reads
    unique_number_aligned_reads = (*aligned_genome_index_data).size();

    // get unique_number_aligned_genomes
    unique_number_aligned_genomes = max_genome_index + 1;

}

//=============================================================================
// Destructor
//=============================================================================
SIGMA_IPOPT_NLP::~SIGMA_IPOPT_NLP()
{
    delete aligned_genome_index_data;
    delete aligned_qvalue_data;
    delete genome_index_info;
	delete reads_info;
	delete reads_weights;
}


//=============================================================================
// Returns the size of the problem 
//=============================================================================
bool SIGMA_IPOPT_NLP::get_nlp_info(Index& n,    // number of variables
        Index& m,           // number of constraints
        Index& nnz_jac_g,   // number of non-zero entries in Jacobian
        Index& nnz_h_lag,   // number of non-zero entries in Hessian
        IndexStyleEnum& index_style)    // C style or Fortran style
{
//echo		std::cout <<"get_nlp_info" << std::endl;

    // number of variables that is the number of "unique_number_aligned_genomes"
    // x[0] through x[n-1]
    n = unique_number_aligned_genomes;

    // one equality constraint (x[0]+x[1]+...+x[n-1] = 1)
    m = 1;

    // in this example the jacobian is dense and contains n nonzeros
    nnz_jac_g = n;

    // the hessian is also dense and has n*n total nonzeros, but we
    // only need the lower left corner (since it is symmetric)
    nnz_h_lag = n*(n+1)/2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}


//=============================================================================
// Returns the variable bounds
//=============================================================================
bool SIGMA_IPOPT_NLP::get_bounds_info(
        Index n,        // number of variables
        Number* x_l,    // lower bound of x
        Number* x_u,    // upper bound of x
        Index m,        // number of constraints
        Number* g_l,    // lower bound of constraints
        Number* g_u)    // upper count of constraints
{
//echo		std::cout <<"get_bounds_info" << std::endl;

    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == unique_number_aligned_genomes);
    assert(m == 1);

    // the variables have lower bounds of 0
    for (Index i=0; i<unique_number_aligned_genomes; i++) {
        x_l[i] = 0.0;
    }

    // the variables have upper bounds of 1
    for (Index i=0; i<unique_number_aligned_genomes; i++) {
        x_u[i] = 1.0;
    }

    // the first constraint g1 has an equality constraint, so we set the
    // upper and lower bound to the same value "1"
    g_l[0] = g_u[0] = 1;

    return true;
}


//=============================================================================
// Returns the initial point for the problem
//=============================================================================
bool SIGMA_IPOPT_NLP::get_starting_point(
        Index n,        // number of variables
        bool init_x,    // (in) if true, must provide an initial value for x
        Number* x,      // (out) the intial value for the primal varialbles, x
        bool init_z,    // (in), if true, provide an inital value for the bound multipliers
        Number* z_L,    // no use
        Number* z_U,    // no use
        Index m,        // (in) the number of constraints
        bool init_lambda,   // (in) true, provide an initial value for the constraint multipliers
        Number* lambda) // no use
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // generate initial point as close as possible to the results
    // map for counting of the highest q-value reads per genome
    std::map<unsigned int, double> count_max_reads_map;

    // for-loop aligned_qvalue_data row
	double total_weights = 0.0;
    for (unsigned int i=0; i<unique_number_aligned_reads; i++ ) {
        unsigned int column_number = (*aligned_qvalue_data)[i].size();
		double read_weight = (*reads_weights)[i];    //modified
		total_weights += read_weight;    //modified

        // initialize max q-value and of the genome index
        double max_qvalue = 0.0;
        unsigned int max_qvalue_genome_index = (*aligned_genome_index_data)[i][0];

        // get max qvalue for each column
        for (unsigned int j=0; j<column_number; j++ ) {
            double this_qvalue = (*aligned_qvalue_data)[i][j];
            unsigned int this_qvalue_genome_index = (*aligned_genome_index_data)[i][j];
            if (this_qvalue > max_qvalue) {
                max_qvalue = this_qvalue;
                max_qvalue_genome_index = this_qvalue_genome_index;
            }
        }

        // count add to count_max_reads_map
        // if not exist
        if (count_max_reads_map.find(max_qvalue_genome_index) == count_max_reads_map.end() ) {
            count_max_reads_map.insert ( std::pair<unsigned int, double>(max_qvalue_genome_index, read_weight) );    //modified
        }
        // else (exist)
        else {
            count_max_reads_map[max_qvalue_genome_index ] += read_weight;    //modified
        }
    }

    // calculate initial x[] by the count of max reads
    std::map<unsigned int, unsigned int>::iterator map_element;
    for (Index j=0; j < unique_number_aligned_genomes; j++ ) {    
        Index aligned_genome_index = j;
        // double percentage_genome = 0.0;
        double probability_genome = 0.0;
        if (count_max_reads_map.find(aligned_genome_index) != count_max_reads_map.end() ) {
            probability_genome = ((double)count_max_reads_map[aligned_genome_index])/((double)total_weights);    //modified
        }
        x[aligned_genome_index] = probability_genome;
        //std::cout << "*** x_0["  << aligned_genome_index << "] = " << x[aligned_genome_index] << std::endl;
    }

    return true;
}


//=============================================================================
// Returns the value of the objective function
//=============================================================================
bool SIGMA_IPOPT_NLP::eval_f(
        Index n,            // number of variables
        const Number* x,    // (in), the values for the primal variables
        bool new_x,         // (in), false if any evaluation method was previous called with the same values in x
        Number& obj_value)  // (out), the value of the objective function
{
//echo	std::cout <<"eval_f" << std::endl;
    // IpOpt internal variable for object function
    obj_value = 0.; 

    // initialize
    unsigned int genome_index;
    double qvalue;

    // objective function:
    // minimize { - sum   ( log( sum   (Q[i][j]*x[j]) ) ) }
    //              i=1..l       j=1..n

    // for loop reads
    int i, column_number, j;
    double obj_value_one;
    #pragma omp parallel for schedule(dynamic) \
    private(i,obj_value_one,column_number,j,genome_index,qvalue)

    for (i=0; i< (int) unique_number_aligned_reads; i++) {

        // save partial objective value for each row
        obj_value_one = 0.; 

        // column number is the number of genomes that matched to the read
        column_number = (*aligned_qvalue_data)[i].size();    
 
        // loop columns of the row
        for (j=0; j<column_number; j++ ) {    

            // get genome index and qvalue from data
            genome_index = (*aligned_genome_index_data)[i][j];
            qvalue = (*aligned_qvalue_data)[i][j];

            // one row of object value inside log is sum{Q[i][j]*X[j]} for j
            obj_value_one += qvalue*x[genome_index]; 
        }   

        // object value = -sum{log[obj_value_one]} for i
        #pragma omp atomic
        obj_value += (-log(obj_value_one)) * (*reads_weights)[i];     //modified
    }
  
    return true;
}


//=============================================================================
// Return the gradient of the objective function grad_{x} f(x)
//=============================================================================
bool SIGMA_IPOPT_NLP::eval_grad_f(
        Index n,            // number of variables
        const Number* x,    // (in), the values for the primal variables
        bool new_x,         // (in), false if any evaluation method was previous called with the same values in x
        Number* grad_f)     // (out), the array of values for the gradient of the objective function (Delta f(x))
{
//echo	std::cout <<"eval_grad_f" << std::endl;
    // initailize grad_f[] to 0.0
    for (Index i=0; i<unique_number_aligned_genomes; i++) {
        grad_f[i] = 0.0;
    }

    // initialize variables
    int genome_index, i, column_number, j;
    double qvalue, grad_f_one, grad_f_denominator_one;

    // gradien of the objective function:
    // grad_f = [ df/dx[0], df/dx[1], .., df/dx[n-1] ]
    // where   df    -Q[0][0]
    //       ----- = ---------------------------------------------------- 
    //       dx[0]   Q[0][0]*x[0] + Q[0][1]*x[1] + ... + Q[0][n-1]*x[n-1]
    //
    //               -Q[1][0]
    //             + ---------------------------------------------------- 
    //               Q[1][0]*x[0] + Q[1][1]*x[1] + ... + Q[1][n-1]*x[n-1]
    //            
    //               ...
    //
    //               -Q[l][0]
    //             + ---------------------------------------------------- 
    //               Q[l][0]*x[0] + Q[l][1]*x[1] + ... + Q[l][n-1]*x[n-1]
    //
  
    // loop reads (from 1 to l)
    #pragma omp parallel for schedule(dynamic) \
    private(i,column_number,grad_f_denominator_one,genome_index,qvalue,j,grad_f_one)
    for (i=0; i< (int) unique_number_aligned_reads; i++) {

        // column number is the number of genomes that matched to the read
        column_number = (*aligned_qvalue_data)[i].size();    
        grad_f_denominator_one = 0.;

        // loop columns
        for (j=0; j<column_number; j++ ) {    

            // get genome index and qvalue from data
            genome_index = (*aligned_genome_index_data)[i][j];
            qvalue = (*aligned_qvalue_data)[i][j];

            // add denominator (Q[i][j]*X[j]) for each mapped genome
            grad_f_denominator_one += qvalue*x[genome_index];
        }

        // loop columns again
        for (j=0; j<column_number; j++ ) {    

            // get genome index and qvalue from data
            genome_index = (*aligned_genome_index_data)[i][j];
            qvalue = (*aligned_qvalue_data)[i][j];

            // calculate grad_f_one for each genome, and sum to grad_f
            grad_f_one = (double)(-qvalue / grad_f_denominator_one) * (*reads_weights)[i];    //modified
			//std::cout << (double)(-qvalue / grad_f_denominator_one) << std::endl;
            #pragma omp atomic
            grad_f[genome_index] += grad_f_one;
        }
    }

    return true;
}


// return the value of the constraints: g(x)
bool SIGMA_IPOPT_NLP::eval_g(Index n, // (in), number of variables
        const Number* x, // (in), the values for the primal variables
        bool new_x, // (in), false if any evaluation method was previous called with the same values in x
        Index m, // (in), the number of constraints
        Number* g) // (out), the array of constraint function values, g(x)
{
//echo		std::cout <<"eval_g" << std::endl;

    // initailize g[0] = 0.0
    // we only have one constraint, so only g[0]
    g[0] = 0.0;

    // Constraints:
    // g[0] = x[0] + x[1] + ... + x[n-1]

    // loop genomes
    for (Index j=0; j<unique_number_aligned_genomes; j++) {
        g[0] += x[j];
    }

    return true;
}


//=============================================================================
// Return the structure or values of the jacobian
//=============================================================================
bool SIGMA_IPOPT_NLP::eval_jac_g(
        Index n,            // (in), number of variables
        const Number* x,    // (in), the values for the primal variables
        bool new_x,         // (in), false if any evaluation method was previous called with the same values in x
        Index m,            // (in), the number of constraints
        Index nele_jac,     // (in), the number of nonzero element
        Index* iRow,        // (out), the row indices of entries in the Jacobian of the constraints
        Index *jCol,        // (out), the column indices of entries in the Jacobian of the constraints
        Number* values)     // (out), the values of the entries in the Jacobian of the constraints
{
//echo		std::cout <<"eval_jac_g" << std::endl;

    if (values == NULL) {
        // return the structure of the jacobian
        for (Index j=0; j<unique_number_aligned_genomes; j++) {
            iRow[j] = 0;
            jCol[j] = j;
        }
    }
    // jac_g[0] = 1
    // jac_g[1] = 1
    // ...
    // jac_g[n-1] = 1
    else {
        // return the values of the jacobian of the constraints
        for (Index j=0; j<unique_number_aligned_genomes; j++) {
            values[j] = 1.0;
        }
    }

    return true;
}


//=============================================================================
// Return the structure or values of the hessian
//=============================================================================
bool SIGMA_IPOPT_NLP::eval_h(
        Index n,                // (in), number of variables
        const Number* x,        // (in), the values for the primal variables
        bool new_x,             // (in), false if any evaluation method was previous called with the same values in x
        Number obj_factor,      // (in), factor in front of the objective term in the Hessian, (sigma)
        Index m,                // (in), the number of constraints
        const Number* lambda,   // (in), the values for the constraint multipliers
        bool new_lambda,        // (in), false if any evaluation method was previous called with the same values in lambda
        Index nele_hess,        // (in), the number of non-zero elements in the Hessian
        Index* iRow,            // (out), the row indices of entries in the Hessian
        Index* jCol,            // (out), the column indices of entries in the Hessian
        Number* values)         // (out), the values of the entries in the Hessian
{
//echo		std::cout <<"eval_h" << std::endl;

    // initialize row and column indices of entries in the Hessian
    if (values == NULL) {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        // the hessian for this problem is actually dense
        Index idx=0;
        for (Index row = 0; row < unique_number_aligned_genomes; row++) {
            for (Index col = 0; col <= row; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

        assert(idx == nele_hess);
    }
    else {
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only

        // column number is the number of genomes that matched to the reads
        int column_number;
        double hess_f_numerator_one;
        double hess_f_denominator_one;

        unsigned int row_genome_index, col_genome_index;
        unsigned int values_idx;
        int i, j, k;

        // Step1: initialize values
        for (Index row = 0; row < unique_number_aligned_genomes; row++) {
            for (Index col = 0; col <= row; col++) {
                values_idx = (row*(row+1)/2)+col;
                values[values_idx] = 0.0;
            }
        }

        // Step2: loop reads to calculate each column for the read
        #pragma omp parallel for schedule(dynamic) \
        private(i,column_number,hess_f_numerator_one,hess_f_denominator_one,j, \
                row_genome_index,col_genome_index,k,values_idx)

        for (i=0; i< (int) unique_number_aligned_reads; i++) {

            // column number is the number of genomes that matched to the reads
            column_number = (*aligned_qvalue_data)[i].size();    
            // initialize for each read
            hess_f_numerator_one = 0.;
            hess_f_denominator_one = 0.;

            // loop columns
            for (j=0; j<column_number; j++ ) {    

                // get genome index and qvalue from data
                //genome_index = (*aligned_genome_index_data)[i][j];
                //qvalue = (*aligned_qvalue_data)[i][j];

                // calculate Q[i][j]*X[j]
                hess_f_denominator_one += (*aligned_qvalue_data)[i][j]*x[(*aligned_genome_index_data)[i][j]];
            }
            hess_f_denominator_one = hess_f_denominator_one*hess_f_denominator_one;


            // loop columns for row of Hessian
            for (j=0; j<column_number; j++ ) {    

                // get genome index and qvalue from data
                row_genome_index = (*aligned_genome_index_data)[i][j];

                // loop columns for columns of Hessian
                for (k=0; k<column_number; k++ ) {    

                    // get genome index and qvalue from data
                    col_genome_index = (*aligned_genome_index_data)[i][k];

                    // only consider low triangle of Hessian
                    if (row_genome_index >= col_genome_index) {

                        // calculate numerator Q[i][j] * Q[p][q]
                        hess_f_numerator_one = (*aligned_qvalue_data)[i][j]*(*aligned_qvalue_data)[i][k];

                        // calculate values_index 
                        values_idx = (row_genome_index*(row_genome_index+1)/2)+col_genome_index;

                        // calculate partial hessian and sum up
                        #pragma omp atomic
                        values[values_idx] += (hess_f_numerator_one / hess_f_denominator_one)*(*reads_weights)[i];     //modified
                    }
                }
            } 
        } 

        // Step3: Multiply obj factor
        for (Index row = 0; row < unique_number_aligned_genomes; row++) {
            for (Index col = 0; col <= row; col++) {
                Index values_idx = (row*(row+1)/2)+col;
                values[values_idx] = (obj_factor * values[values_idx]);  
            }
        }
    }

    return true;
}

//=============================================================================
// This method is called by IPOPT after the algorithm has finished
//=============================================================================
void SIGMA_IPOPT_NLP::finalize_solution(
        SolverReturn status,    // (in), gives the status of the algorithm
        Index n,                // (in), number of variables
        const Number* x,        // (in), the values for the primal variables
        const Number* z_L,      // (in), the final values for the lower bound multiplier
        const Number* z_U,      // (in), the final values for the upper bound multiplier
        Index m,                // (in), the number of constraints
        const Number* g,        // (in), the final value of the constraint function values, g(x)
        const Number* lambda,   // (in), the final value of the constraint multipliers, lambda
        Number obj_value,       // (in), the final value of the objective, f(x)
        const IpoptData* ip_data,           // no use
        IpoptCalculatedQuantities* ip_cq)   // no use
{
    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.

    // print Objective_Value
    std::ofstream output_file (output_filename.c_str());
    std::streamsize old_precision = output_file.precision(20);

    // double
    double percentage_chance_genome;

    // print comments
    output_file << "#\t+\tMatrixName\tTotalNumberReads\tMatchedReads\tUnmatchedReads\tReadLimit" << std::endl;
    output_file << "#\t@\tGenomeIndex\tGenomeName\tNumber\tRelativePercentage" << std::endl;
	
    // print genome index info
	output_file << (*reads_info)[0] << std::endl;
    for (unsigned int i=0; i<(*genome_index_info).size(); i++) {
        percentage_chance_genome = 100*x[i];
        
        output_file << (*genome_index_info)[i] << "\t" << percentage_chance_genome << std::endl;
    }

    // close output_file
    output_file.precision(old_precision);
    output_file.close();
}
