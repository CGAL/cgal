// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_MCS_GET_NORMAL_H
#define CGAL_MCS_GET_NORMAL_H

namespace CGAL {

namespace internal {

template <class Traits>
void normalize(typename Traits::Vector_3& v, const Traits& traits)
{
  double norm = std::sqrt(traits.compute_squared_length_3_object()(v));
  v = traits.construct_divided_vector_3_object()(v, norm);
}


template <class Vertex, class Traits>
typename Traits::Vector_3 get_vertex_normal(
  Vertex& v,
  const Traits& traits)
{
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Vertex::Halfedge_around_vertex_const_circulator HV_circulator;
  Vector normal = traits.construct_vector_3_object()(CGAL::NULL_VECTOR);
  HV_circulator he = v.vertex_begin();
  HV_circulator end = he;
  CGAL_For_all(he,end)
  {
    if(!he->is_border())
    {
      const Point& prev = he->prev()->vertex()->point();
      const Point& curr = he->vertex()->point();
      const Point& next = he->next()->vertex()->point();

      Vector p1 = traits.construct_vector_3_object()(curr, next);
      normalize(p1, traits);
      Vector p2 = traits.construct_vector_3_object()(curr, prev);
      normalize(p2, traits);

      double cosine = traits.compute_scalar_product_3_object()(p1, p2);
      if      (cosine < -1.0) cosine = -1.0;
      else if (cosine >  1.0) cosine =  1.0;
      double angle = acos(cosine);

      Vector n = traits.construct_cross_product_vector_3_object()(
        traits.construct_vector_3_object()(curr, next),
        traits.construct_vector_3_object()(curr, prev) );
      normalize(n, traits);
      n = traits.construct_scaled_vector_3_object()(n, angle);

      normal = traits.construct_sum_of_vectors_3_object()(normal, n);
    }
  }
  normalize(normal, traits);
  return normal;
}

} //namespace internal

} //namespace CGAL

#endif // CGAL_MCS_GET_NORMAL_H
