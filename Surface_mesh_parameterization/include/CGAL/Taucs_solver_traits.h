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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_TAUCS_SOLVER_TRAITS
#define CGAL_TAUCS_SOLVER_TRAITS

#include <CGAL/Taucs_matrix.h>
#include <CGAL/Taucs_vector.h>
#include <CGAL/Taucs_fix.h>

#include <cassert>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>

CGAL_BEGIN_NAMESPACE


/// The class Taucs_symmetric_solver_traits
/// is a traits class for solving SYMMETRIC DEFINITE POSITIVE sparse linear systems
/// using TAUCS solvers family.
/// The default solver is the Multifrontal Supernodal Cholesky Factorization.
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

    /// Create a TAUCS sparse linear solver for SYMMETRIC DEFINITE POSITIVE matrices.
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

#ifdef DEBUG_TRACE
       // Turn on TAUCS trace
       std::cerr.flush();
       taucs_logfile("stderr");
#endif

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

private:

    /// Test if a floating point number is (close to) 0.0.
    static inline bool IsZero(NT a)
    {
        return (CGAL_CLIB_STD::fabs(a) < 10.0 * std::numeric_limits<NT>::min());
    }

// Fields
private:
    const char**  m_options;
    const void**  m_arguments;
};


/// The class Taucs_solver_traits
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

    /// Create a TAUCS sparse linear solver for GENERAL (aka unsymmetric) matrices.
    Taucs_solver_traits()
    {
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

#ifdef DEBUG_TRACE
       // Turn on TAUCS trace
       std::cerr.flush();
       taucs_logfile("stderr");
#endif

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

            // create temporary multifile for out-of-core swapping
            char*   matrixfile = tmpnam(NULL);
            taucs_io_handle* oocL = taucs_io_create_multifile(matrixfile);
            matrixfile = NULL;
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

private:

    /// Test if a floating point number is (close to) 0.0.
    static inline bool IsZero(NT a)
    {
        return (CGAL_CLIB_STD::fabs(a) < 10.0 * std::numeric_limits<NT>::min());
    }

    /// Generate a temporary filename. See "man tmpnam".
    /// This is a replacement for the standard tmpnam() function that is deprecated.
    /// Note: this is a modified version of cupsTempFile().
    char *                                    // Out - Filename or
    tmpnam(char *filename)                    // In - Pointer to buffer
    {
        int           fd;                     // File descriptor for temp file
#ifdef WIN32
        char          tmpdir[1024];           // Windows temporary directory
#else
        char          *tmpdir;                // TMPDIR environment var
#endif
        struct timeval curtime;               // Current time
        static char   buf[L_tmpnam] = "";     // Buffer if you pass in NULL

        // See if a filename was specified.
        if (filename == NULL)
            filename = buf;

        // Get temporary directory.
#ifdef WIN32
        if ((tmpdir = getenv("TEMP")) == NULL)
        {
            GetTempPath(sizeof(tmpdir), tmpdir);
        }
#else
        if ((tmpdir = getenv("TMPDIR")) == NULL)
        {
            // Put root temp files in restricted temp directory.
            if (getuid() == 0)
                tmpdir = "/tmp";
            else
                tmpdir = "/var/tmp";
        }
#endif

        // Make the temporary name using the specified directory.
        do
        {
            // Get the current time of day...
            gettimeofday(&curtime, NULL);

            // Format a string using the hex time values...
            snprintf(filename, L_tmpnam - 1, "%s/%08x%05x", tmpdir,
                    (int)curtime.tv_sec, (int)curtime.tv_usec);

            // Open the file in "exclusive" mode, making sure that we don't
            // stomp on an existing file or someone's symlink crack.
#ifdef O_NOFOLLOW
            fd = open(filename, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0600);
#else
            fd = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0600);
#endif
        }
        while (fd < 0);

        // Close the temp file - it'll be reopened later as needed.
        close(fd);

        // Return the temp filename.
       return (filename);
    }

}; // Taucs_solver_traits


CGAL_END_NAMESPACE

#endif // CGAL_TAUCS_SOLVER_TRAITS
