// Copyright (c) 2020-2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Dmitry Anisimov
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H
#define CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H

#if defined(CGAL_USE_OSQP) || defined(DOXYGEN_RUNNING)

// STL includes.
#include <tuple>
#include <vector>
#include <utility>
#include <exception>

// CGAL includes.
#include <CGAL/assertions.h>

// OSQP includes.
#include <osqp/osqp.h>

namespace CGAL {

/*!
  \ingroup PkgSolverInterfaceNLP

  \brief wraps the external OSQP solver.

  This class provides an interface for formulating and solving
  constrained or unconstrained quadratic programs using the \ref thirdpartyOSQP "OSQP"
  library, which must be available on the system.

  \tparam FT
  number type

  \cgalModels `QuadraticProgramTraits`
*/
template<typename FT>
class OSQP_quadratic_program_traits
{
  using Triplet = std::tuple<std::size_t, std::size_t, FT>; // row, col, value

private:
  std::size_t n; // number of variables
  std::size_t m; // number of constraints

  std::vector<Triplet> P_vec, A_vec;
  std::vector<FT> q_vec, l_vec, u_vec;

public:
  /// Default constructor
  OSQP_quadratic_program_traits() : n(0), m(0) { }

  /// Constructor
  /// \param n the number of variables
  OSQP_quadratic_program_traits(const std::size_t n)
    : n(n), m(0)
  {
    // 6*n, just guessing that this is some kind of sparse problem on a nice 2D mesh
    P_vec.reserve(6 * n);
    q_vec.reserve(n);
  }

  /// Constructor
  /// \param n the number of variables
  /// \param m the number of constraints
  OSQP_quadratic_program_traits(const std::size_t n, const std::size_t m)
    : n(n), m(m)
  {
    P_vec.reserve(6 * n);
    q_vec.reserve(n);

    A_vec.reserve(m);
    l_vec.reserve(m);
    u_vec.reserve(m);
  }

public:
  /// Reset the problem, removing all existing entries and setting sizes to `0`.
  void clear()
  {
    n = m = 0;
    P_vec.clear();
    A_vec.clear();
    q_vec.clear();
    l_vec.clear();
    u_vec.clear();
  }

  /// Change the number of variables and the number of constraints of the problem.
  ///
  /// \warning Calling this function also clears all previous entries.
  void resize(const int new_n,
              const int new_m = 0)
  {
    clear();
    n = new_n;
    m = new_m;
  }

  /// \cond SKIP_IN_MANUAL
  void set_P(const std::size_t i, const std::size_t j, const FT value)
  {
    if(j < i)
      return;
    if(j >= n) // no need to test i since j >= i
      n = j+1;

    P_vec.emplace_back(i, j, value);
  }

  void set_q(const std::size_t i, const FT value)
  {
    if(i >= n)
      n = i+1;
    if(i >= q_vec.size())
      q_vec.resize(n);

    q_vec[i] = value;
  }

  void set_r(const FT) {
    // is not used here
  }

  void set_A(const std::size_t i, const std::size_t j, const FT value)
  {
    if(i >= m)
      m = i+1;
    if(j >= n)
      n = j+1;

    A_vec.emplace_back(i, j, value);
  }

  void set_l(const std::size_t i, const FT value)
  {
    if(i >= m)
      m = i+1;
    if(i >= l_vec.size())
      l_vec.resize(m);

    l_vec[i] = value;
  }

  void set_u(const std::size_t i, const FT value)
  {
    if(i >= m)
      m = i+1;
    if(i >= u_vec.size())
      u_vec.resize(m);

    u_vec[i] = value;
  }

  template<typename OutIterator>
  bool solve(OutIterator solution,
             const FT eps_abs = 1e-10,
             const FT eps_rel = 1e-15,
             const bool verbose = false)
  {
    std::cout << "num variables = " << n << std::endl;
    std::cout << "num constraints = " << m << std::endl;

    CGAL_precondition(n >= 1 && m >= 0);

    CGAL_precondition(q_vec.size() == n);
    CGAL_precondition(l_vec.size() == m && l_vec.size() == m);

    const c_int P_nnz = static_cast<c_int>(P_vec.size());
    c_float P_x[P_nnz];
    c_int   P_i[P_nnz];
    c_int   P_p[n + 1];
    set_matrix_from_triplets("P", P_vec, P_x, P_i, P_p);
    std::cout << "P_nnz: " << P_nnz << std::endl;

    const c_int A_nnz = static_cast<c_int>(A_vec.size());
    c_float A_x[A_nnz];
    c_int   A_i[A_nnz];
    c_int   A_p[n + 1];
    set_matrix_from_triplets("A", A_vec, A_x, A_i, A_p);
    std::cout << "A_nnz: " << A_nnz << std::endl;

    const c_int q_size = static_cast<c_int>(q_vec.size());
    const c_int l_size = static_cast<c_int>(l_vec.size());
    const c_int u_size = static_cast<c_int>(u_vec.size());

    c_float q_x[q_size];
    c_float l_x[l_size];
    c_float u_x[u_size];
    set_qlu_data(q_x, l_x, u_x);

    // Problem settings.
    OSQPSettings *settings = (OSQPSettings *) malloc(sizeof(OSQPSettings));
    CGAL_assertion(settings);

    // Structures.
    OSQPWorkspace *work;
    OSQPData *data = (OSQPData *) malloc(sizeof(OSQPData));
    CGAL_assertion(data);

    // Populate data.
    data->n = static_cast<c_int>(n);
    data->m = static_cast<c_int>(m);
    data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
    CGAL_assertion(data->P);

    data->q = q_x;
    data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
    CGAL_assertion(data->A);

    data->l = l_x;
    data->u = u_x;

    // Set solver settings.
    osqp_set_default_settings(settings);
    settings->eps_abs = eps_abs;
    settings->eps_rel = eps_rel;
    settings->verbose = verbose;

    // Set workspace.
    osqp_setup(&work, data, settings);

    // Solve problem.
    int exitflag = -1;
    try
    {
      exitflag = osqp_solve(work);
    }
    catch (std::exception& e)
    {
      std::cerr << "ERROR: OSQP solver has thrown an exception!" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    const bool success = (exitflag == 0);

    // Create solution.
    const c_float *x = work->solution->x;
    for(std::size_t i=0; i<n; ++i)
    {
      const FT value{x[i]};
      *(++solution) = value;
    }

    // Clean workspace.
    osqp_cleanup(work);
    c_free(data->A);
    c_free(data->P);
    c_free(data);
    c_free(settings);

    return success;
  }
  /// \endcond

public:
  // Based on the code in scipy, function coo_tocsr()
  void set_matrix_from_triplets(const std::string /* name */,
                                const std::vector<Triplet>& triplets,
                                c_float *M_x,
                                c_int *M_i,
                                c_int *M_p) const
  {
    const std::size_t nnz = triplets.size();

    // Compute the number of non-zero entries in each column of the sparse matrix A
    std::fill(M_p, M_p + n, 0);
    for(std::size_t k=0; k<nnz; ++k)
    {
      M_p[std::get<1>(triplets[k])]++;
    }

    // Fill M_p
    c_int cumsum = 0;
    for(std::size_t j=0; j<n; ++j)
    {
      const c_int tmp = M_p[j];
      M_p[j] = cumsum;
      cumsum += tmp;
    }
    M_p[n] = nnz;

    // Write Ai, Ax into M_i, M_x
    for(std::size_t k=0; k<nnz; ++k)
    {
      const c_int col = static_cast<c_int>(std::get<1>(triplets[k]));
      const c_int dest = M_p[col];

      M_i[dest] = static_cast<c_int>(std::get<0>(triplets[k]));
      M_x[dest] = c_float(CGAL::to_double(std::get<2>(triplets[k])));

      M_p[col]++;
    }

    c_int last = 0;
    for(std::size_t j=0; j<=n; ++j)
    {
      const c_int tmp = M_p[j];
      M_p[j] = last;
      last = tmp;
    }

    // std::cout << name + "_x: ";
    // for(std::size_t i=0; i<nnz; ++i)
    //   std::cout << M_x[i] << " ";
    // std::cout << std::endl;

    // std::cout << name + "_i: ";
    // for(std::size_t i=0; i<nnz; ++i)
    //   std::cout << M_i[i] << " ";
    // std::cout << std::endl;

    // std::cout << name + "_p: ";
    // for(std::size_t i=0; i<(n+1); ++i)
    //   std::cout << M_p[i] << " ";
    // std::cout << std::endl;
  }

  void set_qlu_data(c_float *q_x,
                    c_float *l_x,
                    c_float *u_x) const
  {
    for(std::size_t i=0; i<n; ++i)
    {
      q_x[i] = c_float(CGAL::to_double(q_vec[i]));
    }

    for(std::size_t i=0; i<m; ++i)
    {
      l_x[i] = c_float(CGAL::to_double(l_vec[i]));
      u_x[i] = c_float(CGAL::to_double(u_vec[i]));
    }

    // std::cout << "q_x: ";
    // for(std::size_t i=0; i<n; ++i)
    //   std::cout << q_x[i] << " ";
    // std::cout << std::endl;

    // std::cout << "l_x: ";
    // for(std::size_t i=0; i<m; ++i)
    //   std::cout << l_x[i] << " ";
    // std::cout << std::endl;

    // std::cout << "u_x: ";
    // for(std::size_t i=0; i<m; ++i)
    //   std::cout << u_x[i] << " ";
    // std::cout << std::endl;
  }
};

} // namespace CGAL

#endif // CGAL_USE_OSQP or DOXYGEN_RUNNING

#endif // CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H
