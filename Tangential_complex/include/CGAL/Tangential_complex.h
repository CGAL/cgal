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


#ifndef TANGENTIAL_COMPLEX_H
#define TANGENTIAL_COMPLEX_H

#include <CGAL/basic.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>

#include <vector>

namespace CGAL {
  
/// The class Tangential_complex represents a tangential complex
template <
  typename Kernel, 
  int Intrinsic_dimension,
  typename Tr = Regular_triangulation<Regular_triangulation_euclidean_traits<
                  CGAL::Epick_d<Dimension_tag<Intrinsic_dimension> > > >
>
class Tangential_complex
{
  typedef typename Kernel::FT                       FT;
  typedef typename Kernel::Point_d                  Point;
  typedef typename Kernel::Vector_d                 Vector;

  typedef Tr                                        Triangulation;
  typedef typename Triangulation::Point             Tr_point;
  typedef typename Triangulation::Bare_point        Tr_bare_point;
  typedef typename Triangulation::Vertex_handle     Tr_vertex_handle;
  
  typedef typename std::vector<Vector>              Tangent_space_base;

  typedef typename std::vector<Point>               Point_container;
  typedef typename std::vector<Triangulation>       Tr_container;
  typedef typename std::vector<Tangent_space_base>  TS_container;

  // Stores the index of the original Point in the ambient space
  struct Tr_point_with_index 
  : public Tr_point
  {
    Tr_point_with_index(const Tr_point &p, std::size_t i)
      : Tr_point(p), index(i) {}

    std::size_t index;
  };

public:
  /// Constructor
  Tangential_complex() {}
  
  /// Constructor for a range of points
  template <typename InputIterator>
  Tangential_complex(InputIterator first, InputIterator last)
  : m_points(first, last) {}

  /// Destructor
  ~Tangential_complex() {}

  void compute_tangential_complex()
  {
    Point_container::const_iterator it_p = m_points.begin();
    Point_container::const_iterator it_p_end = m_points.end();
    // For each point p in ambient space
    for (std::size_t i = 0 ; it_p != it_p_end ; ++it_p, ++i)
    {
      m_triangulations.push_back(Triangulation(Intrinsic_dimension));
      Triangulation &local_tr = m_triangulations.back();

      // Estimate the tangent space
      Tangent_space_base ts = compute_tangential_space(*it_p);
      
      // Insert p
      Tr_point wp = project_point(*it_p, ts);
      Tr_point_with_index tpwi(wp, i);
      Tr_vertex_handle vh = local_tr.insert(tpwi);

      // Build a minimal triangulation in the tangent space
      // (we only need the star of p)
      Point_container::const_iterator it2_p = m_points.begin();
      for (std::size_t j = 0 ; it2_p != it_p_end ; ++it2_p, ++j)
      {
        // ith point = p, which is already inserted
        if (j != i)
        {
          Tr_point wp = project_point(*it2_p, ts);
          Tr_point_with_index tpwi(wp, j);
          local_tr.insert_if_in_star(tpwi, vh);
        }
      }
    }
  }

private:
  Tangent_space_base compute_tangential_space(const Point &p) const
  {
    Tangent_space_base ts;
    ts.reserve(Intrinsic_dimension);
    // CJTODO: this is only for a sphere in R^3
    Vector n = Kernel().point_to_vector_d_object()(p); // CJTODO: change that?
    Vector t1(-p[1] - p[2], p[0], p[0]);
    Vector t2(p[1] * t1[2] - p[2] * t1[1],
              p[2] * t1[0] - p[0] * t1[2],
              p[0] * t1[1] - p[1] * t1[0]);

    Kernel k;
    Get_functor<Kernel, Squared_length_tag>::type sqlen(k);
    //Get_functor<Kernel, Scaled_vector_tag>::type scale(k);
    //ts.push_back(scale(t1, 1./CGAL::sqrt(sqlen(t1))));
    //ts.push_back(scale(t2, 1./CGAL::sqrt(sqlen(t2))));

    FT t1_len = CGAL::sqrt(sqlen(t1));
    FT t2_len = CGAL::sqrt(sqlen(t2));
    for (int i = 0 ; i < Ambient_dimension<Vector>::value ; ++i)
    {
      t1[i] /= t1_len;
      t2[i] /= t2_len;
    }
    ts.push_back(t1);
    ts.push_back(t2);

    return ts;
  }

  Tr_point project_point(const Point &p, const Tangent_space_base &ts) const
  {
    std::vector<FT> coords;
    coords.reserve(Intrinsic_dimension);
    for (std::size_t i = 0 ; i < Intrinsic_dimension ; ++i)
    {
      //coords[i] = Kernel().point_to_vector_d_object()(p) * ts[i]; // CJTODO: use that
      Kernel k;
      Get_functor<Kernel, Scalar_product_tag>::type scp(k);
      coords.push_back(scp(k.point_to_vector_d_object()(p), ts[i]));
    }

    return Tr_point(
      Tr_bare_point(Intrinsic_dimension, coords.begin(), coords.end()), 0); // CJTODO: poids
  }

private:
  Point_container     m_points;
  TS_container        m_tangent_spaces;
  Tr_container        m_triangulations;

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
