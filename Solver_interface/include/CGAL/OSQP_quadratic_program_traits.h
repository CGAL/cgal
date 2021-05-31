// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H
#define CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H

#if defined(CGAL_USE_OSQP) || defined(DOXYGEN_RUNNING)

// STL includes.
#include <tuple>
#include <vector>
#include <utility>

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
  class OSQP_quadratic_program_traits {
    // row, col, value
    using Triplet = std::tuple<std::size_t, std::size_t, FT>;

  public:
    /// \cond SKIP_IN_MANUAL
    void reserve_P(const std::size_t k) {
      P_vec.reserve(k);
    }

    void reserve_q(const std::size_t n) {
      q_vec.reserve(n);
    }

    void reserve_A(const std::size_t k) {
      A_vec.reserve(k);
    }

    void reserve_l(const std::size_t m) {
      l_vec.reserve(m);
    }

    void reserve_u(const std::size_t m) {
      u_vec.reserve(m);
    }

    void set_P(const std::size_t i, const std::size_t j, const FT value) {
      P_vec.push_back(std::make_tuple(i, j, value));
    }

    void set_q(const std::size_t, const FT value) {
      q_vec.push_back(value);
    }

    void set_r(const FT) {
      // is not used here
    }

    void set_A(const std::size_t i, const std::size_t j, const FT value) {
      A_vec.push_back(std::make_tuple(i, j, value));
    }

    void set_l(const std::size_t, const FT value) {
      l_vec.push_back(value);
    }

    void set_u(const std::size_t, const FT value) {
      u_vec.push_back(value);
    }

    std::size_t P_size() const {
      return P_vec.size();
    }

    std::size_t q_size() const {
      return q_vec.size();
    }

    std::size_t A_size() const {
      return A_vec.size();
    }

    std::size_t l_size() const {
      return l_vec.size();
    }

    std::size_t u_size() const {
      return u_vec.size();
    }

    template<typename OutIterator>
    bool solve(OutIterator solution) {

      const std::size_t num_cols = std::get<1>(P_vec.back()) + 1;
      CGAL_precondition(q_vec.size() == num_cols);
      CGAL_precondition(l_vec.size() == u_vec.size());

      const c_int P_nnz = static_cast<c_int>(P_vec.size());
      c_float P_x[P_nnz];
      c_int   P_i[P_nnz];
      c_int   P_p[P_nnz + 1];
      set_matrix_from_triplets("P", P_vec, P_x, P_i, P_p);
      // std::cout << "P_nnz: " << P_nnz << std::endl;

      const c_int A_nnz = static_cast<c_int>(A_vec.size());
      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[P_nnz + 1];
      set_matrix_from_triplets("A", A_vec, A_x, A_i, A_p);
      // std::cout << "A_nnz: " << A_nnz << std::endl;

      const c_int q_size = static_cast<c_int>(q_vec.size());
      const c_int l_size = static_cast<c_int>(l_vec.size());
      const c_int u_size = static_cast<c_int>(u_vec.size());

      c_float q_x[q_size];
      c_float l_x[l_size];
      c_float u_x[u_size];
      set_qlu_data(q_x, l_x, u_x);

      // Problem settings.
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures.
      OSQPWorkspace *work;
      OSQPData *data;

      // Populate data.
      const c_int n = q_size; // number of variables
      const c_int m = l_size; // number of constraints
      // std::cout << "n: " << n << std::endl;
      // std::cout << "m: " << m << std::endl;

      data = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(n, n, P_nnz, P_x, P_i, P_p);
      data->q = q_x;
      data->A = csc_matrix(m, n, A_nnz, A_x, A_i, A_p);
      data->l = l_x;
      data->u = u_x;

      // Set solver settings.
      osqp_set_default_settings(settings);
      settings->eps_rel = 1.0e-15;
      settings->verbose = false;

      // Set workspace.
      osqp_setup(&work, data, settings);

      // Solve problem.
      const int exitflag = osqp_solve(work);
      const bool success = exitflag == 0 ? true : false;

      // Create solution.
      const c_float *x = work->solution->x;
      for (c_int i = 0; i < n; ++i) {
        const FT value = static_cast<FT>(x[i]);
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

  private:
    std::vector<Triplet> P_vec, A_vec;
    std::vector<FT> q_vec, l_vec, u_vec;

    void set_matrix_from_triplets(
      const std::string /* name */,
      const std::vector<Triplet>& triplets,
      c_float *M_x, c_int *M_i, c_int *M_p) const {

      M_p[0] = 0;
      std::size_t count = 0;
      const std::size_t num_cols = std::get<1>(triplets.back());
      std::size_t ref_col = 0;
      for (ref_col = 0; ref_col <= num_cols; ++ref_col) {
        std::size_t num_rows = 0;
        for (std::size_t i = 0; i < triplets.size(); ++i) {
          const std::size_t row = std::get<0>(triplets[i]);
          const std::size_t col = std::get<1>(triplets[i]);
          const double value = CGAL::to_double(std::get<2>(triplets[i]));

          if (col == ref_col) {
            M_i[count] = row; M_x[count] = value;
            ++count; ++num_rows;
          }
        }
        M_p[ref_col + 1] = M_p[ref_col] + num_rows;
      }

      // std::cout << name + "_x: ";
      // for (std::size_t i = 0; i < count; ++i)
      //   std::cout << M_x[i] << " ";
      // std::cout << std::endl;

      // std::cout << name + "_i: ";
      // for (std::size_t i = 0; i < count; ++i)
      //   std::cout << M_i[i] << " ";
      // std::cout << std::endl;

      // std::cout << name + "_p: ";
      // for (std::size_t i = 0; i <= ref_col; ++i)
      //   std::cout << M_p[i] << " ";
      // std::cout << std::endl;
    }

    void set_qlu_data(
      c_float *q_x, c_float *l_x, c_float *u_x) const {

      CGAL_assertion(q_vec.size() > 1);
      CGAL_assertion(l_vec.size() == u_vec.size());
      const std::size_t n = q_vec.size(); // number of variables
      const std::size_t m = l_vec.size(); // number of constraints

      CGAL_assertion(n > 1);
      for (std::size_t i = 0; i < n; ++i)
        q_x[i] = CGAL::to_double(q_vec[i]);

      CGAL_assertion(m >= 0);
      for (std::size_t i = 0; i < m; ++i) {
        l_x[i] = CGAL::to_double(l_vec[i]);
        u_x[i] = CGAL::to_double(u_vec[i]);
      }

      // std::cout << "q_x: ";
      // for (std::size_t i = 0; i < n; ++i)
      //   std::cout << q_x[i] << " ";
      // std::cout << std::endl;

      // std::cout << "l_x: ";
      // for (std::size_t i = 0; i < m; ++i)
      //   std::cout << l_x[i] << " ";
      // std::cout << std::endl;

      // std::cout << "u_x: ";
      // for (std::size_t i = 0; i < m; ++i)
      //   std::cout << u_x[i] << " ";
      // std::cout << std::endl;
    }
  };

} // namespace CGAL

#endif // CGAL_USE_OSQP or DOXYGEN_RUNNING

#endif // CGAL_OSQP_QUADRATIC_PROGRAM_TRAITS_H
