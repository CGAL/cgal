// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Mario Botsch

#ifndef CGAL_TAUCS_SOLVER_H
#define CGAL_TAUCS_SOLVER_H

#include <CGAL/Taucs_fix.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

// Uncomment the next line to see libraries selected by auto-link
//#define CGAL_LIB_DIAGNOSTIC
#include <CGAL/auto_link/TAUCS.h>

#ifdef WIN32
  #include <CGAL/Win32_exception.h>
#endif

#include <vector>
#include <set>


CGAL_BEGIN_NAMESPACE


// Forward declaration
template<class T> struct Taucs_traits;


/// CLASS Taucs_solver:
/// direct solver for symmetric positive definite sparse systems.
template<class T>       // Tested with T = float or double
class Taucs_solver
{
// Public types
public:

  typedef std::vector<T> Vector;
  typedef T NT;


// Public operations
public:

  Taucs_solver(unsigned int size, ///< Number of rows = number of columns
               unsigned int nb_elements_per_line=0) ///< Number of non null elements per line
                                                    ///< (ignored if 0).
    : PAP(0),
      L(0),
      SL(0),
      n_rows(size),
      row_index(0),
      perm(0),
      invperm(0),
      supernodal_(true),
      m_io_handle(NULL)
  {
    // Reserve taucs_ccs_matrix's arrays to avoid memory fragmentation:
    // - rowind[] = array of non null elements of the matrix, ordered by columns.
    // - values[] = array of row index of each element of rowind[].
    // - colptr[j] is the index of the first element of the column j (or where it
    //   should be if it doesn't exist) + the past-the-end index of the last column.
    colptr.reserve(n_rows+1);
    if (nb_elements_per_line > 0)
    {
      rowind.reserve(n_rows*nb_elements_per_line + 16); // 16 is slack
      values.reserve(n_rows*nb_elements_per_line + 16); // 16 is slack
    }

#ifdef DEBUG_TRACE
    // Turn on TAUCS trace to stderr or to a log file
  #if DEBUG_TRACE >= 2
    std::cerr.flush();
    taucs_logfile((char*)"stderr");
  #else
    taucs_logfile((char*)"taucs.log");
  #endif
#endif
  }

  ~Taucs_solver()
  {
    delete_matrices();
  }

  void begin_row()
  {
    if (colptr.empty() || colptr.back() != (int)values.size())
    {
      colptr.push_back((int)values.size());
      row_index = (int)colptr.size()-1;
    }
    CGAL_surface_reconstruction_points_assertion(row_index <= n_rows);
  }

  void add_value(int _i, T _val)
  {
      // We store only the lower diagonal matrix
      if(_i <= row_index)
      {
        // print warning if reallocating the arrays
        if (values.capacity() == values.size() || rowind.capacity() == rowind.size())
          std::cerr << "Taucs_solver: reallocating the matrix\n";

        values.push_back(_val);
        rowind.push_back(_i);
      }
  }

  void end_row()
  {
    if (colptr.empty() || colptr.back() != (int)values.size())
    {
      colptr.push_back((int)values.size());
      row_index = (int)(colptr.size()-1);
    }
    CGAL_surface_reconstruction_points_assertion(row_index <= n_rows + 1);
  }

  bool factorize(bool _use_supernodal = true)
  {
#ifdef WIN32
    Win32_exception_handler eh; // catch Win32 structured exceptions
#endif

    supernodal_ = _use_supernodal;

    // delete old matrices
    delete_matrices();

    CGAL_surface_reconstruction_points_assertion(perm == NULL);
    CGAL_surface_reconstruction_points_assertion(invperm == NULL);
    CGAL_surface_reconstruction_points_assertion(PAP == NULL);
    CGAL_surface_reconstruction_points_assertion(SL == NULL);
    CGAL_surface_reconstruction_points_assertion(L == NULL);

    finalize_matrix();

    // bandlimitation
    try {
      taucs_ccs_order(&A, &perm, &invperm, (char*)"metis");
    }
    catch (...) {}
    if (perm == NULL || invperm == NULL)
    {
      std::cerr << "Taucs_solver: metis failed\n";
      return false;
    }

    try {
      PAP = taucs_ccs_permute_symmetrically(&A, perm, invperm);
    }
    catch (...) {}
    if (PAP == NULL)
    {
      std::cerr << "Taucs_solver: permutation failed\n";
      return false;
    }

    // Cholesky factorization
    try {
      if (supernodal_)  SL = taucs_ccs_factor_llt_mf (PAP);
      else               L = taucs_ccs_factor_llt    (PAP, 0, 0);
    }
    catch (...) {}
    if (!(L || SL))
    {
      std::cerr << "Taucs_solver: factorization failed\n";
      return false;
    }

    return true;
  }

  bool factorize_ooc()
  {
#ifdef WIN32
    Win32_exception_handler eh; // catch Win32 structured exceptions
#endif

    // delete old matrices
    delete_matrices();

    CGAL_surface_reconstruction_points_assertion(perm == NULL);
    CGAL_surface_reconstruction_points_assertion(invperm == NULL);
    CGAL_surface_reconstruction_points_assertion(PAP == NULL);
    CGAL_surface_reconstruction_points_assertion(m_io_handle == NULL);

    finalize_matrix();

    // bandlimitation
    try {
      taucs_ccs_order(&A, &perm, &invperm, (char*)"metis");
    }
    catch (...) {}
    if (perm == NULL || invperm == NULL)
    {
      std::cerr << "Taucs_solver: metis failed\n";
      return false;
    }

    try {
      PAP = taucs_ccs_permute_symmetrically(&A, perm, invperm);
    }
    catch (...) {}
    if (PAP == NULL)
    {
      std::cerr << "Taucs_solver: permutation failed\n";
      return false;
    }

    // out-of-core Cholesky factorization
    char* matrixfile = tempnam(NULL, "taucs.L");
    CGAL_surface_reconstruction_points_assertion(matrixfile != NULL);
    try {
      m_io_handle = taucs_io_create_multifile(matrixfile);
    }
    catch (...) {}
    free(matrixfile);
    if(m_io_handle == NULL)
    {
        std::cerr << "Taucs_solver: ooc file opening failed\n";
        return false;
    }

    int result = TAUCS_ERROR;
    try {
      double available_memory = taucs_available_memory_size();
      result = taucs_ooc_factor_llt(PAP,m_io_handle,available_memory);
    }
    catch (...) {}
    if(result != TAUCS_SUCCESS)
    {
      std::cerr << "Taucs_solver: ooc factorization failed\n";
      return false;
    }

    return true;
  }



  bool solve(Vector& _b, Vector& _x)
  {
    const unsigned int N = A.n;

    if (N != _b.size() || N != _x.size())
    {
      std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
      return false;
    }

    Vector  PB(N), PX(N);

    // permute rhs
    for (unsigned int i=0; i<N; ++i)
      PB[i] = _b[perm[i]];

    // solve by back-substitution
    if ((supernodal_ ?
        taucs_supernodal_solve_llt(SL, &PX[0], &PB[0]) :
        taucs_ccs_solve_llt(L, &PX[0], &PB[0]))
        != TAUCS_SUCCESS)
    {
      std::cerr << "Taucs_solver: back-substitution failed\n";
      return false;
    }

    // re-permute x
    for (unsigned int i=0; i<N; ++i)
      _x[i] = PX[invperm[i]];

    return true;
  }

  bool solve_ooc(Vector& _b, Vector& _x)
  {
    const unsigned int N = A.n;

    if (N != _b.size() || N != _x.size())
    {
      std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
      return false;
    }

    Vector  PB(N), PX(N);

    // permute rhs
    for (unsigned int i=0; i<N; ++i)
      PB[i] = _b[perm[i]];

    if(taucs_ooc_solve_llt(m_io_handle,&PX[0],&PB[0]) != TAUCS_SUCCESS)
    {
      std::cerr << "Taucs_solver: back-substitution failed\n";
      return false;
    }

    // re-permute x
    for (unsigned int i=0; i<N; ++i)
      _x[i] = PX[invperm[i]];

    return true;
  }


  bool solve_minres(Vector& _b,
                    Vector& _x)
  {
    T* x = new T[_x.size()];
    T* b = new T[_b.size()];
    for(unsigned int i=0;i<_b.size();i++)
    {
      x[i] = _x[i];
      b[i] = _b[i];
    }
    int result = taucs_minres(&A,NULL,NULL,x,b,1000,0.01);
    delete [] x;
    delete [] b;
    return (result == TAUCS_SUCCESS);
  }

  /// Solve using a TAUCS sparse linear solver for symmetric positive definite matrices.
  bool solve_linear(Vector& _b,
                    Vector& _x)
  {
    T* x = new T[_x.size()];
    T* b = new T[_b.size()];
    for(unsigned int i=0;i<_b.size();i++)
    {
      x[i] = _x[i];
      b[i] = _b[i];
    }
    // char* options[] = { "taucs.factor.LU=true", NULL };
    char* options[] = {NULL};
    int result = taucs_linsolve(&A,NULL,1,x,b,options,NULL);
    delete [] x;
    delete [] b;
    return (result == TAUCS_SUCCESS);
  }

  bool solve_conjugate_gradient(Vector& b,
                                Vector& x,
                                const int itermax,
                                const T convergetol)
  {
    finalize_matrix();
    int result = taucs_conjugate_gradients(&A,NULL,NULL,&x[0],&b[0],itermax,convergetol);
    return (result == TAUCS_SUCCESS);
  }


  bool solve(Vector& b,
             Vector& x,
             const unsigned int dim)
  {
    const unsigned int N = A.n;
    CGAL_surface_reconstruction_points_assertion(dim >= 1);

    if (N != b.size()/dim || N != x.size()/dim)
    {
      std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
      return false;
    }

    Vector  PB(N);
    Vector  PX(N);

    // solver component-wise
    for (unsigned int c=0;c<dim;c++)
    {
      // permute rhs
      for(unsigned int i=0; i<N; ++i)
        PB[i] = b[c*N+perm[i]];

      // solve by back-substitution
      if((supernodal_ ?
          taucs_supernodal_solve_llt(SL, &PX[0], &PB[0]) :
          taucs_ccs_solve_llt(L, &PX[0], &PB[0])) != TAUCS_SUCCESS)
      {
        std::cerr << "Taucs_solver: back-substitution failed\n";
        return false;
      }

      // re-permute x
      for (unsigned int i=0; i<N; ++i)
        x[c*N + i] = PX[invperm[i]];
    }
    return true;
  }

// Private operations
private:

  /// Copy constructor and operator =() are not implemented.
  Taucs_solver(const Taucs_solver& toCopy);
  Taucs_solver& operator =(const Taucs_solver& toCopy);

  // Setup ccs matrix.
  void finalize_matrix()
  {
    A.n        = (int)(colptr.size()-1);
    A.m        = (int)(colptr.size()-1);

    // Convert matrix's T type to the corresponding TAUCS constant
    A.flags    = (Taucs_traits<T>::TAUCS_FLAG | TAUCS_SYMMETRIC | TAUCS_LOWER);

    A.colptr   = &colptr[0];
    A.rowind   = &rowind[0];
    A.values.v = &values[0];
  }

  void delete_matrices()
  {
    if (PAP)      { taucs_ccs_free(PAP);               PAP     = 0; }
    if (L)        { taucs_ccs_free(L);                 L       = 0; }
    if (SL)       { taucs_supernodal_factor_free(SL);  SL      = 0; }
    if (perm)     { free(perm);                        perm    = 0; }
    if (invperm)  { free(invperm);                     invperm = 0; }

    if(m_io_handle != NULL)
    {
      taucs_io_delete(m_io_handle);
      // taucs_io_close(m_io_handle);
      m_io_handle = NULL;
    }
  }

// Data
private:

  taucs_ccs_matrix           A, *PAP, *L;
  void                       *SL;
  std::vector<T>             values;
  std::vector<int>           colptr;
  std::vector<int>           rowind;
  int                        n_rows; // number of rows = number of columns
  int                        row_index; // index of current row
  int                        *perm, *invperm;
  bool                       supernodal_;
  taucs_io_handle            *m_io_handle;
};


// Utility class:
// convert matrix's T type to the corresponding TAUCS constant (called TAUCS_FLAG).
template<class T> struct Taucs_traits {};
template<> struct Taucs_traits<double> {
    enum { TAUCS_FLAG = TAUCS_DOUBLE };
};
template<> struct Taucs_traits<float>  {
    enum { TAUCS_FLAG = TAUCS_SINGLE };
};
template<> struct Taucs_traits<taucs_dcomplex> {
    enum { TAUCS_FLAG = TAUCS_DCOMPLEX };
};
template<> struct Taucs_traits<taucs_scomplex> {
    enum { TAUCS_FLAG = TAUCS_SCOMPLEX };
};


CGAL_END_NAMESPACE

#endif // CGAL_TAUCS_SOLVER_H
