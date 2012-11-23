// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_TAUCS_SOLVER_TRAITS_H
#define CGAL_TAUCS_SOLVER_TRAITS_H

#include <CGAL/basic.h> // include basic.h before testing #defines

// Uncomment the next line to see libraries selected by auto-link
//#define CGAL_LIB_DIAGNOSTIC
#include <CGAL/auto_link/TAUCS.h>

#include <CGAL/Taucs_matrix.h>
#include <CGAL/Taucs_vector.h>
#include <CGAL/Taucs_fix.h>

#ifdef WIN32
  #include <CGAL/Win32_exception.h>
#endif
    
#include <boost/shared_ptr.hpp>

#include <stdio.h> // For tempnam
#include <cmath>
#include <cfloat>
#include <climits>

#ifdef _MSC_VER
#include <io.h>
#endif

namespace CGAL {

/// @cond SKIP_IN_MANUAL

/// \ingroup  PkgSurfaceParameterizationAlgebra
///
/// The class Taucs_symmetric_solver_traits
/// is a traits class for solving symmetric positive definite sparse linear systems
/// using TAUCS solvers family.
/// The default solver is the Multifrontal Supernodal Cholesky Factorization.
///
/// \cgalModels `SparseLinearAlgebraTraits_d`

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

    /// Create a TAUCS sparse linear solver for symmetric positive definite matrices.
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
    /// \pre A.row_dimension()    == B.dimension().
    /// \pre A.column_dimension() == X.dimension().
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;          // TAUCS does not support homogeneous coordinates

#ifdef DEBUG_TRACE
        // Turn on TAUCS trace to stderr or to a log file
    #if DEBUG_TRACE >= 2
        std::cerr.flush();
        taucs_logfile((char*)"stderr");
    #else
        taucs_logfile((char*)"taucs.log");
    #endif

//         // Print A and B
//         int n = A.row_dimension();
//         if (n < 20)	// if small matrix, print it entirely
//         {
//           fprintf(stderr, "******************  A:  ******************\n");
//           for (int i=0; i<n; i++)  {
//             for (int j=0; j<n; j++)
//               fprintf(stderr, "%lf\t", (double)A.get_coef(i, j));
//             fprintf(stderr, "\n");
//           }
//           fprintf(stderr, "******************  B:  ******************\n");
//           for (int j=0; j<n; j++)
//             fprintf(stderr, "%lf\t", (double)B[j]);
//           fprintf(stderr, "\n");
//           fprintf(stderr, "******************************************\n");
//         }
//         else		// if large matrix, print only not null elements
//         {
//           fprintf(stderr, "******************  A*X=B  ******************\n");
//           for (int i=0; i<n; i++)  {
//             for (int j=0; j<n; j++)
//               if ( ! IsZero(A.get_coef(i, j)) )
//                 fprintf(stderr, "A[%d][%d] = %lf\t", i, j, (double)A.get_coef(i, j));
//             fprintf(stderr, "\n");
//           }
//           for (int j=0; j<n; j++)
//             if ( ! IsZero(B[j]) )
//               fprintf(stderr, "B[%d] = %lf\t", j, (double)B[j]);
//           fprintf(stderr, "\n");
//           fprintf(stderr, "******************************************\n");
//         }
#endif

#ifdef WIN32
        Win32_exception_handler eh; // catch Win32 structured exceptions
#endif
    
        try
        {
//printf("A[0][0]=%lf\n", (double) A.get_coef(0,0));
//printf("A[77][77]=%lf\n", (double) A.get_coef(77,77));
//printf("taucs_linsolve()\n");
            // Factor, solve and free
            int success = taucs_linsolve((taucs_ccs_matrix*) A.get_taucs_matrix(),
                                         NULL,
                                         1,
                                         X.get_taucs_vector(),
                                         (T*) B.get_taucs_vector(),
                                         (char**) m_options,
                                         (void**) m_arguments);
//printf("A[0][0]=%lf\n", (double) A.get_coef(0,0));
//printf("A[77][77]=%lf\n", (double) A.get_coef(77,77));
            if (success != TAUCS_SUCCESS) {
                taucs_printf((char*)"\tSolving Failed\n");
                return false;
            } else {
                return true;
            }
        }
        catch (...)
        {
            taucs_printf((char*)"\tIncorrect Matrix\n");
            return false;
        }
    }

private:

    // Test if a floating point number is (close to) 0.0.
    static inline bool IsZero(NT a)
    {
        return (CGAL::abs(a) < 10.0 * (std::numeric_limits<NT>::min)());
    }

// Fields
private:
    const char**  m_options;
    const void**  m_arguments;
};


/// \ingroup  PkgSurfaceParameterizationAlgebra
///
/// The class Taucs_solver_traits
/// is a traits class for solving general, that is symmetric and unsymmetric, sparse linear systems
/// using TAUCS out-of-core LU factorization.
///
/// \cgalModels `SparseLinearAlgebraTraits_d`

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

    /// Create a TAUCS sparse linear solver for GENERAL (aka unsymmetric) matrices.
    Taucs_solver_traits()
    {
    }

    /// Solve the sparse linear system "A*X = B".
    /// Return true on success. The solution is then (1/D) * X.
    ///
    /// \pre A.row_dimension()    == B.dimension().
    /// \pre A.column_dimension() == X.dimension().
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;          // TAUCS does not support homogeneous coordinates

#ifdef DEBUG_TRACE
        // Turn on TAUCS trace to stderr or to a log file
    #if DEBUG_TRACE >= 2
        std::cerr.flush();
        taucs_logfile((char*)"stderr");
    #else
        taucs_logfile((char*)"taucs.log");
    #endif

//         // Print A and B
//         int n = A.row_dimension();
//         if (n < 20)  // if small matrix, print it entirely
//         {
//           fprintf(stderr, "******************  A:  ******************\n");
//           for (int i=0; i<n; i++)  {
//             for (int j=0; j<n; j++)
//               fprintf(stderr, "%lf\t", (double)A.get_coef(i, j));
//             fprintf(stderr, "\n");
//           }
//           fprintf(stderr, "******************  B:  ******************\n");
//           for (int j=0; j<n; j++)
//             fprintf(stderr, "%lf\t", (double)B[j]);
//           fprintf(stderr, "\n");
//           fprintf(stderr, "******************************************\n");
//         }
//         else     // if large matrix, print only not null elements
//         {
//           fprintf(stderr, "******************  A*X=B  ******************\n");
//           for (int i=0; i<n; i++)  {
//             for (int j=0; j<n; j++)
//               if ( ! IsZero(A.get_coef(i, j)) )
//                 fprintf(stderr, "A[%d][%d] = %lf\t", i, j, (double)A.get_coef(i, j));
//             fprintf(stderr, "\n");
//           }
//           for (int j=0; j<n; j++)
//             if ( ! IsZero(B[j]) )
//               fprintf(stderr, "B[%d] = %lf\t", j, (double)B[j]);
//           fprintf(stderr, "\n");
//           fprintf(stderr, "******************************************\n");
//         }
#endif

#ifdef WIN32
        Win32_exception_handler eh; // catch Win32 structured exceptions
#endif
    
        try
        {
            int     success;

            // ordering
            int*    perm_raw = NULL;
            int*    invperm_raw = NULL;
            taucs_ccs_order((taucs_ccs_matrix*) A.get_taucs_matrix(),
                            &perm_raw,
                            &invperm_raw,
                            (char*)"colamd");
            boost::shared_ptr<int> perm(perm_raw, free);
            boost::shared_ptr<int> invperm(invperm_raw, free);
            if ( perm == NULL || invperm == NULL)
                throw std::runtime_error("Ordering Failed");

            // Create multi-file for out-of-core swapping.
            // Note: g++ complains that tempnam() is deprecated. You may safely ignore the warning.
#ifdef _MSC_VER
            char template_name[13] = {'t', 'a', 'u', 'c', 's','.','X','X','X','X','X','X', '\0' };
            char* matrixfile = _mktemp(template_name);
            if (matrixfile == NULL)
                throw std::runtime_error("Cannot Create Multifile");
            boost::shared_ptr<taucs_io_handle> oocL(taucs_io_create_multifile(matrixfile), taucs_io_delete);
#else
            boost::shared_ptr<char> matrixfile(tempnam(NULL, "taucs.L"), free);
            if (matrixfile == NULL)
                throw std::runtime_error("Cannot Create Multifile");
            boost::shared_ptr<taucs_io_handle> oocL(taucs_io_create_multifile(matrixfile.get()), taucs_io_delete);
#endif
            if (oocL == NULL)
                throw std::runtime_error("Cannot Create Multifile");

            // factor
            int memory_mb = int(taucs_available_memory_size()/1048576.0);
            success = taucs_ooc_factor_lu((taucs_ccs_matrix*) A.get_taucs_matrix(),
                                           perm.get(),
                                           oocL.get(),
                                           memory_mb*1048576.0);
            if (success != TAUCS_SUCCESS)
                throw std::runtime_error("Factorization Failed");

            // solve
            success = taucs_ooc_solve_lu(oocL.get(),
                                         X.get_taucs_vector(),
                                        (T*) B.get_taucs_vector());
            if (success != TAUCS_SUCCESS)
                throw std::runtime_error("Solving Failed");

            return true;
        }
        catch (std::exception& e)
        {
            taucs_printf((char*)"\t");
            taucs_printf((char*)(e.what() != NULL ? e.what() : "Incorrect Matrix"));
            taucs_printf((char*)"\n");
            return false;
        }
        catch (...)
        {
            taucs_printf((char*)"\tIncorrect Matrix\n");
            return false;
        }
    }

private:

    // Test if a floating point number is (close to) 0.0.
    static inline bool IsZero(NT a)
    {
        return ( ::CGAL::abs(a) < 10.0 * (std::numeric_limits<NT>::min)());
    }
};
/// @endcond

} //namespace CGAL

#endif // CGAL_TAUCS_SOLVER_TRAITS_H
