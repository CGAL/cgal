// Copyright (c) 2024  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDING_BOX_FILTER_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDING_BOX_FILTER_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include<CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/property_map.h>

#include <optional>

namespace CGAL {
namespace Surface_mesh_simplification {


template<class BaseFilter = internal::Dummy_filter>
class Bounding_box_filter
{
public:
  Bounding_box_filter(const BaseFilter& base_filter = BaseFilter())
    : m_base_filter(base_filter)
  {}


  template <typename Profile>
  std::optional<typename Profile::Point>
  operator()(const Profile& profile, std::optional<typename Profile::Point> op) const
  {
    typedef typename Profile::VertexPointMap                              Vertex_point_map;

    typedef typename Profile::Geom_traits                                 Geom_traits;
    typedef typename Geom_traits::Vector_3                                Vector;

    typedef typename boost::property_traits<Vertex_point_map>::value_type Point;
    typedef typename boost::property_traits<Vertex_point_map>::reference  Point_reference;

    const Geom_traits& gt = profile.geom_traits();
    const Vertex_point_map& vpm = profile.vertex_point_map();

    op = m_base_filter(profile, op);
    if(op)
    {
      const Point& placement = *op;
      Bbox_3 bb;
      for(auto vd : profile.link()){
        Point_reference p = get(vpm, vd);
        bb += p.bbox();
      }
      double wx = bb.xmax() - bb.xmin();
      double wy = bb.ymax() - bb.ymin();
      double wz = bb.zmax() - bb.zmin();
      bb = Bbox_3(bb.xmin()-  wx/10.0, bb.ymin() - wy/10.0, bb.zmin()- wz/10.0, bb.xmax() + wx/10.0, bb.ymax() + wy/10.0, bb.zmax()+ wz/10.0);
      if(!do_overlap(bb, placement.bbox())){
             return std::optional<typename Profile::Point>();
      }
    }
    return op;
  }



private:
  const BaseFilter m_base_filter;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDING_BOX_FILTER_H
