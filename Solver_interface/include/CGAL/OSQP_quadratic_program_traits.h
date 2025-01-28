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
#include <memory>

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
  number type that `FieldNumberType`

  \note The `FT` type is provided for convenience. Internally, this FT type is converted
  to `OSQPFloat` type that can be set either to `float` or `double`. By default, the `double`
  type is used. After the optimization is complete, the `OSQPFloat` type is converted back to `FT`.
  See more about `OSQPFloat` <a href="https://osqp.org/docs/interfaces/C.html#data-types">here</a>.

  \cgalModels{QuadraticProgramTraits}
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
  /// %Default constructor
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
  /// Resets the problem, removing all existing entries and setting sizes to `0`.
  void clear()
  {
    n = m = 0;
    P_vec.clear();
    A_vec.clear();
    q_vec.clear();
    l_vec.clear();
    u_vec.clear();
  }

  /// Changes the number of variables and the number of constraints of the problem.
  ///
  /// \warning Calling this function also clears all previous entries.
  void resize(const std::size_t new_n,
              const std::size_t new_m = 0)
  {
    clear();
    n = new_n;
    m = new_m;
    P_vec.reserve(6 * n);
    q_vec.reserve(n);
    if(m > 0)
    {
      A_vec.reserve(m);
      l_vec.reserve(m);
      u_vec.reserve(m);
    }
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
             const double eps_abs = 1e-10,
             const double eps_rel = 1e-15,
             const bool verbose = false)
  {
    if(verbose)
    {
      std::cout << "num variables = " << n << std::endl;
      std::cout << "num constraints = " << m << std::endl;
    }
    CGAL_precondition(n >= 1); // m >= 0

    CGAL_precondition(q_vec.size() == n);
    CGAL_precondition(l_vec.size() == m && l_vec.size() == m);

    const OSQPInt P_nnz = static_cast<OSQPInt>(P_vec.size());
    auto P_x = std::make_unique<OSQPFloat[]>(P_nnz);
    auto P_i = std::make_unique<OSQPInt[]>(P_nnz);
    auto P_p = std::make_unique<OSQPInt[]>(n + 1);
    set_matrix_from_triplets("P", P_vec, P_x.get(), P_i.get(), P_p.get());
    if(verbose) std::cout << "P_nnz: " << P_nnz << std::endl;

    const OSQPInt A_nnz = static_cast<OSQPInt>(A_vec.size());
    auto A_x = std::make_unique<OSQPFloat[]>(A_nnz);
    auto A_i = std::make_unique<OSQPInt[]>(A_nnz);
    auto A_p = std::make_unique<OSQPInt[]>(n + 1);
    set_matrix_from_triplets("A", A_vec, A_x.get(), A_i.get(), A_p.get());
    if(verbose) std::cout << "A_nnz: " << A_nnz << std::endl;

    const OSQPInt q_size = static_cast<OSQPInt>(q_vec.size());
    const OSQPInt l_size = static_cast<OSQPInt>(l_vec.size());
    const OSQPInt u_size = static_cast<OSQPInt>(u_vec.size());

    auto q_x = std::make_unique<OSQPFloat[]>(q_size);
    auto l_x = std::make_unique<OSQPFloat[]>(l_size);
    auto u_x = std::make_unique<OSQPFloat[]>(u_size);
    set_qlu_data(q_x.get(), l_x.get(), u_x.get());

    // Problem settings.
    OSQPSettings *settings = (OSQPSettings *) malloc(sizeof(OSQPSettings));
    CGAL_assertion(settings);

    // Structures.
    OSQPCscMatrix* P = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));
    OSQPCscMatrix* A = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));

    csc_set_data(A, m, n, A_nnz, A_x.get(), A_i.get(), A_p.get());
    csc_set_data(P, n, n, P_nnz, P_x.get(), P_i.get(), P_p.get());

    // Set solver settings.
    osqp_set_default_settings(settings);
    settings->eps_abs = eps_abs;
    settings->eps_rel = eps_rel;
    settings->verbose = verbose;

    OSQPSolver* solver = NULL;
    OSQPInt exitflag = osqp_setup(&solver, P, q_x.get(), A, l_x.get(), u_x.get(), m, n, settings);

    try
    {
      if (!exitflag)
        exitflag = osqp_solve(solver);

      const OSQPFloat* x = solver->solution->x;
      for (std::size_t i = 0; i < n; ++i)
      {
        const FT value{ x[i] };
        *(++solution) = value;
      }
    }
    catch (std::exception& e)
    {
      std::cerr << "ERROR: OSQP solver has thrown an exception!" << std::endl;
      std::cerr << e.what() << std::endl;
    }

    osqp_cleanup(solver);
    if (A) free(A);
    if (P) free(P);
    if (settings) free(settings);

    return (exitflag == 0);
  }
  /// \endcond

private:
  // Based on the code in scipy, function coo_tocsr()
  void set_matrix_from_triplets(const std::string /* name */,
                                const std::vector<Triplet>& triplets,
                                OSQPFloat *M_x,
                                OSQPInt *M_i,
                                OSQPInt *M_p) const
  {
    const std::size_t nnz = triplets.size();

    // Compute the number of non-zero entries in each column of the sparse matrix A
    std::fill(M_p, M_p + n, 0);
    for(std::size_t k=0; k<nnz; ++k)
    {
      M_p[std::get<1>(triplets[k])]++;
    }

    // Fill M_p
    OSQPInt cumsum = 0;
    for(std::size_t j=0; j<n; ++j)
    {
      const OSQPInt tmp = M_p[j];
      M_p[j] = cumsum;
      cumsum += tmp;
    }
    M_p[n] = nnz;

    // Write Ai, Ax into M_i, M_x
    for(std::size_t k=0; k<nnz; ++k)
    {
      const OSQPInt col = static_cast<OSQPInt>(std::get<1>(triplets[k]));
      const OSQPInt dest = M_p[col];

      M_i[dest] = static_cast<OSQPInt>(std::get<0>(triplets[k]));
      M_x[dest] = OSQPFloat(CGAL::to_double(std::get<2>(triplets[k])));

      M_p[col]++;
    }

    OSQPInt last = 0;
    for(std::size_t j=0; j<=n; ++j)
    {
      const OSQPInt tmp = M_p[j];
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

  void set_qlu_data(OSQPFloat *q_x,
                    OSQPFloat *l_x,
                    OSQPFloat *u_x) const
  {
    for(std::size_t i=0; i<n; ++i)
    {
      q_x[i] = OSQPFloat(CGAL::to_double(q_vec[i]));
    }

    for(std::size_t i=0; i<m; ++i)
    {
      l_x[i] = OSQPFloat(CGAL::to_double(l_vec[i]));
      u_x[i] = OSQPFloat(CGAL::to_double(u_vec[i]));
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
