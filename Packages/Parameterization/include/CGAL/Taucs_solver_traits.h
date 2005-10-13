// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_TAUCS_SOLVER_TRAITS
#define CGAL_TAUCS_SOLVER_TRAITS

#include <CGAL/Taucs_matrix.h>
#include <CGAL/Taucs_vector.h>

#include <cassert>
#include <stdio.h>

CGAL_BEGIN_NAMESPACE


/// Class Taucs_symmetric_solver_traits
/// is a traits class for solving SYMMETRIC DEFINIE POSITIVE sparse linear systems
/// using TAUCS solvers family.
/// The default solver is the Multifrontal Supernodal Cholesky Factorization.
///
/// @todo Add to Taucs_symmetric_solver_traits the ability to solve
/// linear systems in the least squares sense.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d concept.

template<class T>       // Tested with T = taucs_single or taucs_double
                        // May also work with T = taucs_dcomplex and taucs_scomplex
class Taucs_symmetric_solver_traits
{
// Public types
public:

    typedef Taucs_symmetric_matrix<T>   Matrix;
    typedef Taucs_vector<T>             Vector;
    typedef T                           NT;

// Public operations
public:

    /// Create a TAUCS sparse linear solver for SYMMETRIC DEFINIE POSITIVE matrices.
    /// The default solver is the Multifrontal Supernodal Cholesky Factorization.
    /// See taucs_linsolve() documentation for the meaning of the
    /// 'options' and 'arguments' parameters.
    Taucs_symmetric_solver_traits(
                    const char*  options[]   = NULL,  ///< must be persistent
		    const void*  arguments[] = NULL)  ///< must be persistent
    {
        static const char* MULTIFRONTAL_LLT[] = {"taucs.factor.LLT=true",
                                                 "taucs.factor.mf=true",
                                                 "taucs.factor.ordering=metis",
                                                 NULL};
        m_options   = (options == NULL) ? MULTIFRONTAL_LLT : options;
        m_arguments = arguments;
    }

    /// Solve the sparse linear system "A*X = B".
    /// Return true on success. The solution is then (1/D) * X.
    ///
    /// Preconditions:
    /// - A.row_dimension()    == B.dimension().
    /// - A.column_dimension() == X.dimension().
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;          // TAUCS does not support homogeneous coordinates

//#ifndef NDEBUG
//        // Turn on TAUCS trace
//        std::cerr.flush();
//        taucs_logfile("stderr");
//#endif

//#ifndef NDEBUG
//        // Debug trace
//        fprintf(stderr, "\n");
//        fprintf(stderr, "linear_solver:\n");
//        int n = A.row_dimension();
//        if (n < 20)	// if small matrix, print it entirely
//        {
//	    fprintf(stderr, "******************  A:  ******************\n");
//	    for (int i=0; i<n; i++)  {
//		    for (int j=0; j<n; j++)
//			    fprintf(stderr, "%lf\t", (double)A.get_coef(i,j));
//		    fprintf(stderr, "\n");
//	    }
//	    fprintf(stderr, "******************  B:  ******************\n");
//	    for (int j=0; j<n; j++)
//		    fprintf(stderr, "%lf\t", (double)B[j]);
//	    fprintf(stderr, "\n");
//	    fprintf(stderr, "******************************************\n");
//        }
//        else		// if large matrix, print only not null elements
//        {
//	    fprintf(stderr, "******************  A*X=B  ******************\n");
//	    for (int i=0; i<n; i++)  {
//		for (int j=0; j<n; j++)
//		    if ( ! IsZero(A.get_coef(i,j)) )
//			fprintf(stderr, "A[%d][%d] = %lf\t", i, j, (double)A.get_coef(i,j));
//		fprintf(stderr, "\n");
//	    }
//	    for (int j=0; j<n; j++)
//		if ( ! IsZero(B[j]) )
//		    fprintf(stderr, "B[%d] = %lf\t", j, (double)B[j]);
//	    fprintf(stderr, "\n");
//	    fprintf(stderr, "******************************************\n");
//        }
//#endif

        try
        {
            // Factor, solve and free
            int success = taucs_linsolve((taucs_ccs_matrix*) A.get_taucs_matrix(),
                                        NULL,
                                        1,
                                        X.get_taucs_vector(),
                                        (T*) B.get_taucs_vector(),
                                        (char**) m_options,
                                        (void**) m_arguments);
            return (success == TAUCS_SUCCESS);
        }
        catch (...)
        {
            // if incorrect matrix
            return false;
        }
    }

    /// Indicate if the linear system can be solved and if the matrix conditioning is good.
    ///
    /// Preconditions:
    /// - A.row_dimension() == B.dimension().
    ///
    /// @todo Implement Taucs_symmetric_solver_traits::is_solvable() by solving the system,
    /// then checking that | ||A*X||/||B|| - 1 | < epsilon
    bool is_solvable (const Matrix& A, const Vector& B)
    {
        // This feature is not implemented in TAUCS => we do only basic checking
        if (A.row_dimension() != B.dimension())
            return false;

        return true;
    }

private:

    /// Test if a floating point number is (close to) 0.0
    static inline bool IsZero(NT a)
    {
        return (CGAL_CLIB_STD::fabs(a) < 10.0 * std::numeric_limits<NT>::min());
    }

// Fields
private:
    const char**  m_options;
    const void**  m_arguments;
};


/// Class Taucs_solver_traits
/// is a traits class for solving GENERAL (aka unsymmetric) sparse linear systems
/// using TAUCS out-of-core LU factorization.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d concept.

template<class T>       // Tested with T = taucs_single or taucs_double
                        // May also work with T = taucs_dcomplex and taucs_scomplex
class Taucs_solver_traits
{
// Public types
public:

    typedef Taucs_matrix<T>             Matrix;
    typedef Taucs_vector<T>             Vector;
    typedef T                           NT;

// Public operations
public:

    /// Solve the sparse linear system "A*X = B".
    /// Return true on success. The solution is then (1/D) * X.
    ///
    /// Preconditions:
    /// - A.row_dimension()    == B.dimension().
    /// - A.column_dimension() == X.dimension().
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;          // TAUCS does not support homogeneous coordinates

//#ifndef NDEBUG
//        // Turn on TAUCS trace
//        std::cerr.flush();
//        taucs_logfile("stderr");
//#endif

//#ifndef NDEBUG
//        // Debug trace
//        fprintf(stderr, "\n");
//        fprintf(stderr, "linear_solver:\n");
//        int n = A.row_dimension();
//        if (n < 20)	// if small matrix, print it entirely
//        {
//	    fprintf(stderr, "******************  A:  ******************\n");
//	    for (int i=0; i<n; i++)  {
//		    for (int j=0; j<n; j++)
//			    fprintf(stderr, "%lf\t", (double)A.get_coef(i,j));
//		    fprintf(stderr, "\n");
//	    }
//	    fprintf(stderr, "******************  B:  ******************\n");
//	    for (int j=0; j<n; j++)
//		    fprintf(stderr, "%lf\t", (double)B[j]);
//	    fprintf(stderr, "\n");
//	    fprintf(stderr, "******************************************\n");
//        }
//        else		// if large matrix, print only not null elements
//        {
//	    fprintf(stderr, "******************  A*X=B  ******************\n");
//	    for (int i=0; i<n; i++)  {
//		for (int j=0; j<n; j++)
//		    if ( ! IsZero(A.get_coef(i,j)) )
//			fprintf(stderr, "A[%d][%d] = %lf\t", i, j, (double)A.get_coef(i,j));
//		fprintf(stderr, "\n");
//	    }
//	    for (int j=0; j<n; j++)
//		if ( ! IsZero(B[j]) )
//		    fprintf(stderr, "B[%d] = %lf\t", j, (double)B[j]);
//	    fprintf(stderr, "\n");
//	    fprintf(stderr, "******************************************\n");
//        }
//#endif

        try
        {
            int     success;

            // ordering
            int*    perm;
            int*    invperm;
            taucs_ccs_order((taucs_ccs_matrix*) A.get_taucs_matrix(),
                            &perm,
                            &invperm,
                            "colamd");
            if (perm == NULL) {
                taucs_printf("\tOrdering Failed\n");
                return false;
            }

            // create multifile for out-of-core swapping
            char*   matrixfile = tempnam(NULL, "taucs.L");
            taucs_io_handle* oocL = taucs_io_create_multifile(matrixfile);
            free(matrixfile); matrixfile = NULL;
            if (oocL == NULL) {
                taucs_printf("\tCannot Create Multifile\n");
                return false;
            }

            // factor
            int memory_mb = int(taucs_available_memory_size()/1048576.0);
            success = taucs_ooc_factor_lu((taucs_ccs_matrix*) A.get_taucs_matrix(),
                                        perm,
                                        oocL,
                                        memory_mb*1048576.0);
            if (success != TAUCS_SUCCESS) {
                taucs_printf("\tFactorization Failed\n");
                return false;
            }

            // solve
            success = taucs_ooc_solve_lu(oocL,
                                        X.get_taucs_vector(),
                                        (T*) B.get_taucs_vector());
            if (success != TAUCS_SUCCESS) {
                taucs_printf("\tSolving Failed\n");
                return false;
            }

            // free
            taucs_io_delete(oocL);

            return true;
        }
        catch (...)
        {
            // if incorrect matrix
            return false;
        }
    }

    /// Indicate if the linear system can be solved and if the matrix conditioning is good.
    ///
    /// Preconditions:
    /// - A.row_dimension() == B.dimension().
    ///
    /// @todo Implement Taucs_solver_traits::is_solvable() by solving the system,
    /// then checking that | ||A*X||/||B|| - 1 | < epsilon
    bool is_solvable (const Matrix& A, const Vector& B)
    {
        // This feature is not implemented in TAUCS => we do only basic checking
        if (A.row_dimension() != B.dimension())
            return false;

        return true;
    }

private:

    /// Test if a floating point number is (close to) 0.0
    static inline bool IsZero(NT a)
    {
        return (CGAL_CLIB_STD::fabs(a) < 10.0 * std::numeric_limits<NT>::min());
    }
};


CGAL_END_NAMESPACE

#endif // CGAL_TAUCS_SOLVER_TRAITS
