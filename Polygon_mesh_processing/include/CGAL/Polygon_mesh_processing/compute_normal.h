// Copyright (c) 2013 INRIA Sophia-Anitpolis (France).
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
// Author(s)     : Pierre Alliez


#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

namespace CGAL{

namespace Polygon_mesh_processing{

/// \ingroup PkgPolygonMeshProcessing
/// computes the outward unit vector normal to facet `f`.
/// `%Kernel::%FT` should be a model of `FieldWithSqrt`
/// \todo make this function working with a `FaceGraph`
template <class Kernel, class Facet>
typename Kernel::Vector_3 compute_facet_normal(const Facet& f)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Facet::Halfedge_around_facet_const_circulator HF_circulator;
  Vector normal = CGAL::NULL_VECTOR;
  HF_circulator he = f.facet_begin();
  HF_circulator end = he;
  CGAL_For_all(he,end)
  {
    const Point& prev = he->prev()->vertex()->point();
    const Point& curr = he->vertex()->point();
    const Point& next = he->next()->vertex()->point();
    Vector n = CGAL::cross_product(next-curr,prev-curr);
    normal = normal + n;
  }
  return normal / std::sqrt(normal * normal);
}

/// \ingroup PkgPolygonMeshProcessing
/// computes the unit normal at vertex `v` as the average of the normals of incident facets.
/// `%Kernel::%FT` should be a model of `FieldWithSqrt`
/// \todo make this function working with a `FaceGraph`
template <class Kernel, class Vertex>
typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
{
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Vertex::Halfedge_around_vertex_const_circulator HV_circulator;
  Vector normal = CGAL::NULL_VECTOR;
  HV_circulator he = v.vertex_begin();
  HV_circulator end = he;
  CGAL_For_all(he,end)
  {
    if(!he->is_border())
    {
      Vector n = compute_facet_normal<Kernel>(*he->facet());
      normal = normal + (n / std::sqrt(n*n));
    }
  }
  return normal / std::sqrt(normal * normal);
}

} } // end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
