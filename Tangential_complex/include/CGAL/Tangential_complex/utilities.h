// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_TC_UTILITIES_H
#define CGAL_TC_UTILITIES_H

#include <CGAL/basic.h>
#include <CGAL/Dimension.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <set>
#include <vector>
#include <atomic> // CJTODO: this is C++11 => use boost.Atomic (but it's too recent) 
                  // or tbb::atomic (works for doubles, but not officially)

// choose exact integral type for QP solver
#ifdef CGAL_USE_GMP
# include <CGAL/Gmpzf.h>
  typedef CGAL::Gmpzf ET;
#else
# include <CGAL/MP_Float.h>
  typedef CGAL::MP_Float ET;
#endif
//#define CGAL_QP_NO_ASSERTIONS // CJTODO: NECESSARY? http://doc.cgal.org/latest/QP_solver/group__PkgQPSolverFunctions.html#ga1fefbd0436aca0e281f88e8e6cd8eb74



namespace CGAL {
namespace Tangential_complex_ {

  // Provides copy constructors to std::atomic so that
  // it can be used in a vector
  template <typename T>
  struct Atomic_wrapper 
    : public std::atomic<T>
  {
    typedef std::atomic<T> Base;

    Atomic_wrapper() {}
    Atomic_wrapper(const T &t) : Base(t) {}
    Atomic_wrapper(const std::atomic<T> &a) : Base(a.load()) {}
    Atomic_wrapper(const Atomic_wrapper &other) : Base(other.load())
    {}
    
    Atomic_wrapper &operator=(const T &other)
    {
      Base::store(other);
      return *this;
    }
    Atomic_wrapper &operator=(const std::atomic<T> &other)
    {
      Base::store(other.load());
      return *this;
    }
    Atomic_wrapper &operator=(const Atomic_wrapper &other)
    {
      Base::store(other.load());
      return *this;
    }
  };

  /*template <typename T>
  struct Atomic_wrapper
  {
    std::atomic<T> _a;

    Atomic_wrapper()
      :_a()
    {}

    Atomic_wrapper(const std::atomic<T> &other)
      :_a(other.load())
    {}

    Atomic_wrapper(const Atomic_wrapper &other)
      :_a(other._a.load())
    {}

    Atomic_wrapper(const T &other)
      :_a(other)
    {}

    Atomic_wrapper &operator=(const std::atomic<T> &other)
    {
      _a.store(other._a.load());
      return *this;
    }

    Atomic_wrapper &operator=(const Atomic_wrapper &other)
    {
      _a.store(other._a.load());
      return *this;
    }

    Atomic_wrapper &operator=(const T &other)
    {
      _a.store(other);
      return *this;
    }

    operator T() const
    {
      return _a.load();
    }

    operator std::atomic<T>() const
    {
      return _a;
    }
  };*/

  // Modifies v in-place
  template <typename K>
  typename K::Vector_d &
  normalize_vector(
    typename K::Vector_d &v,
    K const& k)
  {
    v = k.scaled_vector_d_object()(
      v, FT(1)/CGAL::sqrt(k.squared_length_d_object()(v)));
    return v;
  }

  template <typename K>
  std::vector<typename K::Vector_d>
  compute_gram_schmidt_basis(
    std::vector<typename K::Vector_d> const& input_basis,
    K const& k)
  {
    typedef typename K::FT            FT;
    typedef typename K::Vector_d      Vector;
    typedef std::vector<Vector>       Basis;

    // Kernel functors
    typename K::Squared_length_d        sqlen      = k.squared_length_d_object();
    typename K::Scaled_vector_d         scaled_vec = k.scaled_vector_d_object();
    typename K::Scalar_product_d        inner_pdct = k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec   = k.difference_of_vectors_d_object();

    Basis output_basis;

    typename Basis::const_iterator inb_it = input_basis.begin();
    typename Basis::const_iterator inb_it_end = input_basis.end();
    for (int i = 0 ; inb_it != inb_it_end ; ++inb_it, ++i)
    {
      Vector u = *inb_it;

      typename Basis::iterator outb_it = output_basis.begin();
      for (int j = 0 ; j < i ; ++j)
      {
        Vector const& ej = *outb_it;
        Vector u_proj = scaled_vec(ej, inner_pdct(u, ej));
        u = diff_vec(u, u_proj);
      }

      output_basis.push_back(normalize_vector(u, k));
    }

    return output_basis;
  }
  
  // CJTODO: use CGAL::Combination_enumerator<int> (cf. Tangential_complex.h)
  // Compute all the k-combinations of elements
  // Output_iterator::value_type must be std::set<std::size_t> >
  template <typename Elements_container, typename Output_iterator>
  void combinations(const Elements_container elements, int k,
                    Output_iterator combinations)
  {
    std::size_t n = elements.size();
    std::vector<bool> booleans(n, false);
    std::fill(booleans.begin() + n - k, booleans.end(), true);
    do
    {
      std::set<std::size_t> combination;
      Elements_container::const_iterator it_elt = elements.begin();
      for (std::size_t i = 0 ; i < n ; ++i, ++it_elt)
      {
        if (booleans[i])
          combination.insert(*it_elt);
      }
      *combinations++ = combination;

    } while (std::next_permutation(booleans.begin(), booleans.end()));
  }

  // P: dual face in Delaunay triangulation (p0, p1,… pn)
  // Q: vertices which are common neighbors of all vertices of P
  template <typename K, typename Point, typename Point_range, typename Vector_range>
  bool does_voronoi_face_and_alpha_tangent_subspace_intersect(
    Point const& center_pt,
    Point_range  const& P,
    Point_range  const& Q,
    Vector_range const& orthogonal_subspace_basis,
    typename K::FT alpha,
    K const& k)
  {
    // Notations:
    // Fv: Voronoi k-face
    // Fd: dual, (D-k)-face of Delaunay (p0, p1,… pn)

    typedef typename K::FT                      FT;
    typedef typename K::Point_d                 Point;
    typedef typename K::Vector_d                Vector;
    
    typename K::Scalar_product_d scalar_pdct = k.scalar_product_d_object();
    typename K::Point_to_vector_d pt_to_vec  = k.point_to_vector_d_object();

    const int ambient_dim = k.point_dimension_d_object()(center_pt);

    std::size_t card_P = P.size();
    std::size_t card_Q = Q.size();
    std::size_t card_OSB = orthogonal_subspace_basis.size();
    std::size_t num_couples_among_P = card_P*(card_P-1)/2;
    std::size_t num_equations = 
      2*num_couples_among_P + card_P*card_Q + 2*card_OSB;
    
    // Linear solver
    typedef CGAL::Quadratic_program<FT> Linear_program;
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;

    Linear_program lp(CGAL::SMALLER, false);
    int current_row = 0;

    //=========== First set of equations ===========
    // For point pi in P
    //   2(p0 - pi).x = p0² - pi² 
    Point const& p0 = center_pt;
    FT p0_dot_p0 = scalar_pdct(pt_to_vec(p0), pt_to_vec(p0));
    for (Point_range::const_iterator it_p = P.begin(), it_p_end = P.end() ; 
         it_p != it_p_end ; ++it_p)
    {
      Point const& pi = *it_p;

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(p0[k] - pi[k]));

      lp.set_b(current_row,
               p0_dot_p0 - scalar_pdct(pt_to_vec(pi), pt_to_vec(pi)));
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    // CJTODO: this code might be useful for Option 1
    /*CGAL::Combination_enumerator<int> pi_pj(2, 0, static_cast<int>(card_P));
    for ( ; !pi_pj.finished() ; ++pi_pj)
    {
      Point const& pi = P[pi_pj[0]];
      Point const& pj = P[pi_pj[1]];
      
      for (int k = 0 ; k < ambient_dim ; ++k)
      {
        FT a = 2*(pi[k] + pj[k]);
        lp.set_a(k, current_row    , -a);
        lp.set_a(k, current_row + 1,  a);
      }

      FT b = scalar_pdct(pi, pi) - scalar_pdct(pj, pj);
      lp.set_b(current_row    , -b);
      lp.set_b(current_row + 1,  b);

      current_row += 2;
    }*/

    //=========== Second set of equations ===========
    // For each point qi in Q
    //  2(qi - p0).x <= qi² - p0²
    for (Point_range::const_iterator it_q = Q.begin(), it_q_end = Q.end() ; 
         it_q != it_q_end ; ++it_q)
    {
      Point const& qi = *it_q;

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(qi[k] - p0[k]));

      lp.set_b(current_row, 
               scalar_pdct(pt_to_vec(qi), pt_to_vec(qi)) - p0_dot_p0);

      ++current_row;
    }

    //=========== Third set of equations ===========
    // For each vector of OSB
    //     bi.x <=  bi.p + alpha
    //    -bi.x <= -bi.p + alpha
    for (Vector_range::const_iterator it_osb = 
           orthogonal_subspace_basis.begin(), 
         it_osb_end = orthogonal_subspace_basis.end() ; 
         it_osb != it_osb_end ; ++it_osb)
    {
      Vector const& bi = *it_osb;

      for (int k = 0 ; k < ambient_dim ; ++k)
      {
        lp.set_a(k, current_row    ,  bi[k]);
        lp.set_a(k, current_row + 1, -bi[k]);
      }

      FT bi_dot_p = scalar_pdct(bi, pt_to_vec(center_pt));
      lp.set_b(current_row    ,  bi_dot_p + alpha);
      lp.set_b(current_row + 1, -bi_dot_p + alpha);

      current_row += 2;
    }

    //=========== Other LP parameters ===========
    lp.set_c(0, 1); // Minimize x[0]

    //=========== Solve =========================
    LP_solution solution = CGAL::solve_linear_program(lp, ET());
    if (solution.solves_linear_program(lp))
      std::cout << solution;
    else
      std::cout << "ERROR\n";

    return solution.solves_linear_program(lp);
  }

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_UTILITIES_H
