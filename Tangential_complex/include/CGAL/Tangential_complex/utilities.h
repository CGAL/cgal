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
// (Gmpzf is not thread-safe)
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
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
  typename K::Vector_d& normalize_vector(typename K::Vector_d& v,
                                         K const& k)
  {
    v = k.scaled_vector_d_object()(
      v, typename K::FT(1)/CGAL::sqrt(k.squared_length_d_object()(v)));
    return v;
  }

  template<typename Kernel>
  struct Basis
  {
    typedef typename Kernel::FT                             FT;
    typedef typename Kernel::Point_d                        Point;
    typedef typename Kernel::Vector_d                       Vector;
    typedef typename std::vector<Vector>::const_iterator    const_iterator;

    Point m_origin; //fixme should probably be (const?) ref ?
    std::vector<Vector> m_vectors;

    Point& origin() { return m_origin; }
    const Point& origin() const { return m_origin; }
    const_iterator begin() const { return m_vectors.begin(); }
    const_iterator end() const { return m_vectors.end(); }
    std::size_t size() const { return m_vectors.size(); }

    Vector& operator[](const std::size_t i) { return m_vectors[i]; }
    const Vector& operator[](const std::size_t i) const { return m_vectors[i]; }
    void push_back(const Vector& v) { m_vectors.push_back(v); }
    void reserve(const std::size_t s) { m_vectors.reserve(s); }

    Basis() { }
    Basis(const Point& p) : m_origin(p) { }
    Basis(const Point& p, const std::vector<Vector>& vectors)
      : m_origin(p), m_vectors(vectors)
    { }
    
#ifdef CGAL_ALPHA_TC
    // Thickening vectors...

    struct Thickening_vector
    {
      Thickening_vector() : alpha_minus(FT(0)), alpha_plus(FT(0)) {}

      Vector vec;
      FT alpha_minus;
      FT alpha_plus;
    };
    typedef std::vector<Thickening_vector> Thickening_vectors;

    Thickening_vectors m_thickening_vectors;


    std::size_t num_thickening_vectors() const 
    {
      return m_thickening_vectors.size();
    }
    Thickening_vectors const& thickening_vectors() const 
    {
      return m_thickening_vectors;
    }
    void add_thickening_vector(Vector const& vec, FT alpha_minus, FT alpha_plus)
    {
      m_thickening_vectors.push_back(
        Thickening_vector(vec, alpha_minus, alpha_plus));
    }

    int dimension() const
    {
      return static_cast<int>(m_vectors.size() + m_thickening_vectors.size());
    }
#else
    int dimension() const
    {
      return static_cast<int>(m_vectors.size());
    }
#endif
  };

  template <typename K>
  Basis<K> compute_gram_schmidt_basis(Basis<K> const& input_basis, K const& k)
  {
    typedef Basis<K>                  Basis;
    typedef typename K::Vector_d      Vector;

    // Kernel functors
    typename K::Scaled_vector_d         scaled_vec = k.scaled_vector_d_object();
    typename K::Scalar_product_d        inner_pdct = k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec   = k.difference_of_vectors_d_object();

    Basis output_basis(input_basis.origin());

    typename Basis::const_iterator inb_it = input_basis.begin();
    typename Basis::const_iterator inb_it_end = input_basis.end();
    for (int i = 0 ; inb_it != inb_it_end ; ++inb_it, ++i)
    {
      Vector u = *inb_it;

      for (int j = 0 ; j < i ; ++j)
      {
        Vector const& ej = input_basis[j];
        Vector u_proj = scaled_vec(ej, inner_pdct(u, ej) / inner_pdct(ej, ej));
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
      typename Elements_container::const_iterator it_elt = elements.begin();
      for (std::size_t i = 0 ; i < n ; ++i, ++it_elt)
      {
        if (booleans[i])
          combination.insert(*it_elt);
      }
      *combinations++ = combination;

    } while (std::next_permutation(booleans.begin(), booleans.end()));
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  template <typename K, typename Point_range, typename Weight_range,
            typename Indexed_point_range, typename Indexed_point_range_2,
            typename Basis>
  bool does_voronoi_face_and_alpha_tangent_subspace_intersect(
      Point_range const& all_points,
      Weight_range const& all_weights,
      std::size_t center_pt_index,
      Indexed_point_range const& P,
      Indexed_point_range_2 const& Q,
      Basis const& orthogonal_subspace_basis,
      typename K::FT alpha,
      K const& k)
  {
    // Notations:
    // Fv: Voronoi k-face
    // Fd: dual, (D-k)-face of Delaunay (p0, p1, ..., pn)

    typedef typename K::FT                      FT;
    typedef typename K::Point_d                 Point;
    typedef typename K::Vector_d                Vector;

    typename K::Scalar_product_d scalar_pdct = k.scalar_product_d_object();
    typename K::Point_to_vector_d pt_to_vec  = k.point_to_vector_d_object();
    typename K::Compute_coordinate_d coord = k.compute_coordinate_d_object();

    Point const& center_pt = all_points[center_pt_index];
    int const ambient_dim = k.point_dimension_d_object()(center_pt);

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
    // For points pi in P
    //   2(p0 - pi).x = p0^2 - w0 - pi^2 + wi
    Point const& p0 = center_pt;
    FT const w0 = all_weights[center_pt_index];
    FT p0_dot_p0 = scalar_pdct(pt_to_vec(p0), pt_to_vec(p0));

    for (typename Indexed_point_range::const_iterator it_p = P.begin(),
                                                      it_p_end = P.end() ;
         it_p != it_p_end ; ++it_p)
    {
      Point const& pi = all_points[*it_p];
      FT const wi = all_weights[*it_p];

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(p0, k) - coord(pi, k)));

      FT pi_dot_pi = scalar_pdct(pt_to_vec(pi), pt_to_vec(pi));
      lp.set_b(current_row, p0_dot_p0 - pi_dot_pi - w0 + wi);
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    // CJTODO: this code might be useful for Option 1
    /*CGAL::Combination_enumerator<int> pi_pj(2, 0, static_cast<int>(card_P));
    for ( ; !pi_pj.finished() ; ++pi_pj)
    {
      Point const& pi = P[pi_pj[0]];
      FT wi = all_weights[pi_pj[0]];
      Point const& pj = P[pi_pj[1]];
      FT wj = all_weights[pi_pj[1]];

      for (int k = 0 ; k < ambient_dim ; ++k)
      {
        FT a = 2*(coord(pi, k) + coord(pj, k));
        lp.set_a(k, current_row    , -a);
        lp.set_a(k, current_row + 1,  a);
      }

      FT b = scalar_pdct(pi, pi) - wi - scalar_pdct(pj, pj) + wj;
      lp.set_b(current_row    , -b);
      lp.set_b(current_row + 1,  b);

      current_row += 2;
    }*/

    //=========== Second set of equations ===========
    // For each point qi in Q
    //  2(qi - p0).x <= qi^2 - wi - p0^2 + w0
    for (typename Indexed_point_range_2::const_iterator it_q = Q.begin(),
                                                        it_q_end = Q.end() ;
         it_q != it_q_end ; ++it_q)
    {
      Point const& qi = all_points[*it_q];
      FT const wi = all_weights[*it_q];

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(qi, k) - coord(p0, k)));

      FT qi_dot_qi = scalar_pdct(pt_to_vec(qi), pt_to_vec(qi));
      lp.set_b(current_row, qi_dot_qi - wi - p0_dot_p0 + w0);

      ++current_row;
    }

    //=========== Third set of equations ===========
    // For each vector bi of OSB, (x-p).bi <= alpha and >= -alpha
    // p is the origin of the basis
    //     bi.x <=  bi.p + alpha
    //    -bi.x <= -bi.p + alpha
    for (typename Basis::const_iterator it_osb =
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

      Point const& basis_origin = orthogonal_subspace_basis.origin();
      FT bi_dot_p = scalar_pdct(bi, pt_to_vec(basis_origin));
      lp.set_b(current_row    ,  bi_dot_p + alpha);
      lp.set_b(current_row + 1, -bi_dot_p + alpha);

      current_row += 2;
    }

    //=========== Other LP parameters ===========
    lp.set_c(0, 1); // Minimize x[0]

    //=========== Solve =========================
    LP_solution solution = CGAL::solve_linear_program(lp, ET());
    bool ret = (solution.status() == CGAL::QP_OPTIMAL);

    return ret;
  }

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_UTILITIES_H
