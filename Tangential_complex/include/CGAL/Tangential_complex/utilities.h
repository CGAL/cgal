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

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <set>
#include <vector>
#include <array>
#include <atomic> // CJTODO: this is C++11 => use boost.Atomic (but it's too recent)
                  // or tbb::atomic (works for doubles, but not officially)


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

    std::size_t m_origin;
    std::vector<Vector> m_vectors;

    std::size_t origin() const { return m_origin; }
    void set_origin(std::size_t o) { m_origin = o; }
    const_iterator begin() const { return m_vectors.begin(); }
    const_iterator end() const { return m_vectors.end(); }
    std::size_t size() const { return m_vectors.size(); }

    Vector& operator[](const std::size_t i) 
    {
      return m_vectors[i]; 
    }
    const Vector& operator[](const std::size_t i) const 
    {
      return m_vectors[i];
    }
    void push_back(const Vector& v) 
    {
      m_vectors.push_back(v);
    }
    void reserve(const std::size_t s)
    {
      m_vectors.reserve(s);
    }

    Basis() { }
    Basis(std::size_t origin) : m_origin(origin) { }
    Basis(std::size_t origin, const std::vector<Vector>& vectors)
      : m_origin(origin), m_vectors(vectors)
    { }
    
#ifdef CGAL_ALPHA_TC
    // Thickening vectors...

    struct Thickening_vector
    {
      Thickening_vector() 
        : alpha_minus(FT(0)), alpha_plus(FT(0)) {}
      Thickening_vector(Vector const& v, FT am, FT ap) 
        : vec(v), alpha_minus(am), alpha_plus(ap) {}

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
    void add_thickening_vector(
      Vector const& vec, FT alpha_minus = FT(0), FT alpha_plus = FT(0))
    {
      m_thickening_vectors.push_back(
        Thickening_vector(vec, alpha_minus, alpha_plus));
    }

    void update_alpha_values_of_thickening_vectors(
      Vector const& vec, Kernel const& k)
    {
      typename Kernel::Scalar_product_d k_scalar_pdct =
        k.scalar_product_d_object();

      for (Thickening_vectors::iterator it_v = m_thickening_vectors.begin(),
                                        it_v_end = m_thickening_vectors.end() ;
          it_v != it_v_end ; ++it_v)
      {
        const FT MARGIN_RATIO = 1.001; // CJTODO TEMP
        FT alpha_i = k_scalar_pdct(it_v->vec, vec);
        if (alpha_i * MARGIN_RATIO > it_v->alpha_plus)
        {
#ifdef CGAL_TC_VERY_VERBOSE
          std::cerr << "OLD alpha+ = " << it_v->alpha_plus << std::endl;
#endif
          it_v->alpha_plus = alpha_i * MARGIN_RATIO;
#ifdef CGAL_TC_VERY_VERBOSE
          std::cerr << "NEW alpha+ = " << it_v->alpha_plus << std::endl;
          std::cerr << "NOT MODIFIED alpha- = " << it_v->alpha_minus << std::endl;
#endif
        }
        else if (alpha_i * MARGIN_RATIO < it_v->alpha_minus)
        {
#ifdef CGAL_TC_VERY_VERBOSE
          std::cerr << "OLD alpha- = " << it_v->alpha_minus << std::endl;
#endif
          it_v->alpha_minus = alpha_i * MARGIN_RATIO;
#ifdef CGAL_TC_VERY_VERBOSE
          std::cerr << "NEW alpha- = " << it_v->alpha_minus << std::endl;
          std::cerr << "NOT MODIFIED alpha+ = " << it_v->alpha_plus << std::endl;
#endif
        }
      }
    }

    FT alpha_minus(std::size_t i) const
    {
      return m_thickening_vectors[i].alpha_minus;
    }
    FT alpha_plus(std::size_t i) const
    {
      return m_thickening_vectors[i].alpha_plus;
    }
    
    // Returns 0 if no thickening vectors
    FT max_absolute_alpha() const
    {
      FT max = FT(0);
      
      for (Thickening_vectors::const_iterator 
             it_v = m_thickening_vectors.begin(),
             it_v_end = m_thickening_vectors.end() ;
           it_v != it_v_end ; 
           ++it_v)
      {
        if (it_v->alpha_plus > max)
          max = it_v->alpha_plus;
        if (-it_v->alpha_minus > max)
          max = -it_v->alpha_minus;
      }

      return max;
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
  
  // Using Gram-Schmidt
  // * If the resulting vector after G-S algorithm is below "sqlen_threshold",
  //   the vector considered linearly dependend to the existing vectors 
  //   and is not added to the basis
  // * Returns true if the vector was added to the basis
  template <typename K>
  bool add_vector_to_orthonormal_basis(
    Basis<K> & basis, typename K::Vector_d const& v, K const& k, 
    typename K::FT sqlen_threshold = typename K::FT(1e-13)
#ifdef CGAL_ALPHA_TC
    , bool add_to_thickening_vectors = false
#endif
    )
  {
    typedef Basis<K>                  Basis;
    typedef typename K::FT            FT;
    typedef typename K::Vector_d      Vector;

    // Kernel functors
    typename K::Scaled_vector_d         scaled_vec = k.scaled_vector_d_object();
    typename K::Scalar_product_d        scalar_pdct = k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec   = k.difference_of_vectors_d_object();

    Vector u = v;
    for (int j = 0 ; j < basis.size() ; ++j)
    {
      Vector const& ej = basis[j];
      Vector u_proj = scaled_vec(ej, scalar_pdct(u, ej) / scalar_pdct(ej, ej));
      u = diff_vec(u, u_proj);
    }
    for (int j = 0 ; j < basis.num_thickening_vectors() ; ++j)
    {
      Vector const& ej = basis.thickening_vectors()[j].vec;
      Vector u_proj = scaled_vec(ej, scalar_pdct(u, ej) / scalar_pdct(ej, ej));
      u = diff_vec(u, u_proj);
    }
    
    FT sqlen_new_v = k.squared_length_d_object()(u);
    bool add_it = (sqlen_new_v > sqlen_threshold);
    if (add_it)
    {
      Vector new_v = scaled_vec(u, FT(1)/CGAL::sqrt(sqlen_new_v));

      // If new_v is small, run the Gram-Schmidt once more to 
      // re-orthogonalize it
      if (sqlen_new_v < 0.01)
      {
        for (int j = 0 ; j < basis.size() ; ++j)
        {
          Vector const& ej = basis[j];
          Vector new_v_proj = scaled_vec(
            ej, scalar_pdct(new_v, ej) / scalar_pdct(ej, ej));
          new_v = diff_vec(new_v, new_v_proj);
        }
        for (int j = 0 ; j < basis.num_thickening_vectors() ; ++j)
        {
          Vector const& ej = basis.thickening_vectors()[j].vec;
          Vector new_v_proj = scaled_vec(
            ej, scalar_pdct(new_v, ej) / scalar_pdct(ej, ej));
          new_v = diff_vec(new_v, new_v_proj);
        }
        sqlen_new_v = k.squared_length_d_object()(new_v);
        new_v = scaled_vec(new_v, FT(1)/CGAL::sqrt(sqlen_new_v));
      }

#ifdef CGAL_ALPHA_TC
      if (add_to_thickening_vectors)
        basis.add_thickening_vector(new_v);
      else
#endif
        basis.push_back(new_v);
    }
    return add_it;
  }

  template<
    typename Kernel, typename Tangent_space_basis,
    typename OutputIteratorPoints, typename OutputIteratorTS>
    bool load_points_from_file(
    const std::string &filename,
    OutputIteratorPoints points,
    OutputIteratorTS tangent_spaces,
    std::size_t only_first_n_points = std::numeric_limits<std::size_t>::max())
  {
    typedef typename Kernel::Point_d    Point;
    typedef typename Kernel::Vector_d   Vector;

    std::ifstream in(filename);
    if (!in.is_open())
    {
      std::cerr << "Could not open '" << filename << "'" << std::endl;
      return false;
    }

    bool contains_tangent_spaces =
      (filename.substr(filename.size() - 3) == "pwt");

    Kernel k;
    Point p;
    int num_ppints;
    in >> num_ppints;

    std::size_t i = 0;
    while (i < only_first_n_points && in >> p)
    {
      *points++ = p;
      if (contains_tangent_spaces)
      {
        Tangent_space_basis tsb(i);
        for (int d = 0 ; d < 2 ; ++d) // CJTODO : pas toujours "2"
        {
          Vector v;
          in >> v;
          tsb.push_back(CGAL::Tangential_complex_::normalize_vector(v, k));
        }
        *tangent_spaces++ = tsb;
      }
      ++i;
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

    return true;
  }

  // 1st line: number of points
  // Then one point per line
  template <typename Kernel, typename Point_range>
  std::ostream &export_point_set(
    Kernel const& k,
    Point_range const& points,
    std::ostream & os,
    const char *coord_separator = " ")
  {
    // Kernel functors
    typename Kernel::Construct_cartesian_const_iterator_d ccci =
      k.construct_cartesian_const_iterator_d_object();

    os << points.size() << "\n";

    typename Point_range::const_iterator it_p = points.begin();
    typename Point_range::const_iterator it_p_end = points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      for (auto it = ccci(*it_p) ; it != ccci(*it_p, 0) ; ++it) // CJTODO: C++11
        os << CGAL::to_double(*it) << coord_separator;

      os << "\n";
    }

    return os;
  }

  template <typename K>
  Basis<K> compute_gram_schmidt_basis(Basis<K> const& input_basis, K const& k)
  {
    typedef Basis<K>                  Basis;

    Basis output_basis(input_basis.origin());

    // Add vector one by one
    typename Basis::const_iterator inb_it = input_basis.begin();
    typename Basis::const_iterator inb_it_end = input_basis.end();
    for (int i = 0 ; inb_it != inb_it_end ; ++inb_it, ++i)
      add_vector_to_orthonormal_basis(output_basis, *inb_it, k);

    return output_basis;
  }

  // Functor to compute stereographic projection from S^3(sphere_radius) to R^3
  template <typename K>
  class Orthogonal_projection
  {
  public:
    typedef typename K::FT      FT;
    typedef typename K::Point_d Point;

    // center_of_projection will be sent to infinity by the projection
    Orthogonal_projection(
      std::array<int, 3> const& selected_coords, K const& k)
      : m_selected_coords(selected_coords), m_k(k) 
    {}

    Point operator()(Point const& p) const
    {
      typedef K::FT         FT;
      typedef K::Point_d    Point;

      typename K::Construct_point_d constr_pt = m_k.construct_point_d_object();
      typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

      std::vector<FT> projected_p;
      projected_p.reserve(3);

      for (int i = 0 ; i < 3 ; ++i)
        projected_p.push_back(coord(p, m_selected_coords[i]));

      return constr_pt(3, projected_p.begin(), projected_p.end());
    }

  private:
    std::array<int, 3> m_selected_coords; // CJTODO C++11
    K const& m_k;
  };

  // Functor to compute radial projection from R^4 to S^3(sphere_radius)
  // The returned point coordinates are expressed with the origin 
  // at `center_of_sphere`, i.e. the new points are all on the sphere 
  // S^3{(0,0,0,0), sphere_radius}
  template <typename K>
  class R4_to_S3_radial_projection
  {
  public:
    typedef typename K::FT      FT;
    typedef typename K::Point_d Point;

    // center_of_projection will be sent to infinity by the projection
    R4_to_S3_radial_projection(
      FT sphere_radius, Point const& center_of_sphere, K const& k)
      : m_sphere_radius(sphere_radius), m_center_of_sphere(center_of_sphere),
      m_k(k) {}

    Point operator()(Point const& p) const
    {
      CGAL_assertion(m_k.point_dimension_d_object()(p) == 4);

      typedef K::FT         FT;
      typedef K::Point_d    Point;
      typedef K::Vector_d   Vector;

      typename K::Translated_point_d transl = m_k.translated_point_d_object();
      typename K::Point_to_vector_d pt_to_vec = m_k.point_to_vector_d_object();
      typename K::Vector_to_point_d vec_to_pt = m_k.vector_to_point_d_object();
      typename K::Squared_length_d sqlen = m_k.squared_length_d_object();
      typename K::Scaled_vector_d scaled_vec = m_k.scaled_vector_d_object();

      Point transl_p = transl(p, scaled_vec(pt_to_vec(m_center_of_sphere), FT(-1)));
      Vector v = pt_to_vec(transl_p);
      v = scaled_vec(v, m_sphere_radius / CGAL::sqrt(sqlen(v)));

      return vec_to_pt(v);
    }

  private:
    FT m_sphere_radius;
    Point m_center_of_sphere;
    K const& m_k;
  };

  // Functor to compute stereographic projection from S^3(sphere_radius) to R^3
  template <typename K>
  class S3_to_R3_stereographic_projection
  {
  public:
    typedef typename K::FT      FT;
    typedef typename K::Point_d Point;

    // center_of_projection will be sent to infinity by the projection
    S3_to_R3_stereographic_projection(
      FT sphere_radius, Point const& center_of_projection, K const& k)
    : m_sphere_radius(sphere_radius), m_center_of_proj(center_of_projection), 
      m_k(k) {}

    Point operator()(Point const& p) const
    {
      CGAL_assertion(m_k.point_dimension_d_object()(p) == 4);

      typedef K::FT         FT;
      typedef K::Point_d    Point;

      typename K::Construct_point_d constr_pt = m_k.construct_point_d_object();
      typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

      std::vector<FT> stereo_proj;
      stereo_proj.reserve(3);

      FT t = (2 * m_sphere_radius + coord(m_center_of_proj, 3))
        / (m_sphere_radius - (coord(p, 3) - coord(m_center_of_proj, 3)));
      for (int i = 0 ; i < 3 ; ++i)
        stereo_proj.push_back(coord(m_center_of_proj, i) + t*(coord(p, i) - coord(m_center_of_proj, i)));

      return constr_pt(3, stereo_proj.begin(), stereo_proj.end());
    }

  private:
    FT m_sphere_radius;
    Point m_center_of_proj;
    K const& m_k;
  };

  // Functor to project R^4 points to R^3
  template <typename K>
  class R4_to_R3_using_radial_then_stereographic_projection
  {
  public:
    typedef typename K::FT       FT;
    typedef typename K::Point_d  Point;

    // sphere_radius and center_of_projection are for the stereographic
    // projection
    R4_to_R3_using_radial_then_stereographic_projection(
      FT sphere_radius, Point const& center_of_sphere,
      Point const& center_of_projection, K const& k)
    : m_R4toS3(sphere_radius, center_of_sphere, k),
      m_S3toR3(sphere_radius, center_of_projection, k),
      m_k(k) {}

    Point operator()(Point const& p) const
    {
      return m_S3toR3(m_R4toS3(p));
    }

  private:
    R4_to_S3_radial_projection<K>         m_R4toS3;
    S3_to_R3_stereographic_projection<K>  m_S3toR3;
    K const&                              m_k;
  };

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

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_UTILITIES_H
