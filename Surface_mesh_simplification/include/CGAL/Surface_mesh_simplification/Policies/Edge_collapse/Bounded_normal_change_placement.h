// Copyright (c) 2017  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/property_map.h>

#include <boost/optional.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class GetPlacement>
class Bounded_normal_change_placement
{
public:
  Bounded_normal_change_placement(const GetPlacement& get_placement = GetPlacement())
    : m_get_placement(get_placement)
  {}

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile) const
  {
    typedef typename Profile::VertexPointMap                              Vertex_point_map;

    typedef typename Profile::Geom_traits                                 Geom_traits;
    typedef typename Geom_traits::Vector_3                                Vector;

    typedef typename boost::property_traits<Vertex_point_map>::value_type Point;
    typedef typename boost::property_traits<Vertex_point_map>::reference  Point_reference;

    const Geom_traits& gt = profile.geom_traits();
    const Vertex_point_map& vpm = profile.vertex_point_map();

    boost::optional<typename Profile::Point> op = m_get_placement(profile);
    if(op)
    {
      // triangles returns the triangles of the star of the vertices of the edge to collapse
      // First the two trianges incident to the edge, then the other triangles
      // The second vertex of each triangle is the vertex that gets placed
       const typename Profile::Triangle_vector& triangles = profile.triangles();
       if(triangles.size() > 2)
       {
         typename Profile::Triangle_vector::const_iterator it = triangles.begin();

         if(profile.left_face_exists())
           ++it;
         if(profile.right_face_exists())
           ++it;

         while(it!= triangles.end())
         {
           const typename Profile::Triangle& t = *it;
           Point_reference p = get(vpm, t.v0);
           Point_reference q = get(vpm, t.v1);
           Point_reference r = get(vpm, t.v2);
           const Point& q2 = *op;

           Vector eqp = gt.construct_vector_3_object()(q, p);
           Vector eqr = gt.construct_vector_3_object()(q, r);
           Vector eq2p = gt.construct_vector_3_object()(q2, p);
           Vector eq2r = gt.construct_vector_3_object()(q2, r);

           Vector n1 = gt.construct_cross_product_vector_3_object()(eqp, eqr);
           Vector n2 = gt.construct_cross_product_vector_3_object()(eq2p, eq2r);

           if(!is_positive(gt.compute_scalar_product_3_object()(n1, n2)))
             return boost::optional<typename Profile::Point>();

           ++it;
         }
       }
    }

    return op;
  }

private:
  const GetPlacement m_get_placement;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H
