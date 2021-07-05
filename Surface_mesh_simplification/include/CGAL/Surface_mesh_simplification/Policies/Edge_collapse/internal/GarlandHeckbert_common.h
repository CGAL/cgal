// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_INTERNAL_COMMON_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_INTERNAL_COMMON_H

#include <CGAL/license/Surface_mesh_simplification.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

namespace common {
  template<typename GeomTraits, typename VPM, typename TM>
  typename GeomTraits::Vector_3 construct_unit_normal_from_face(
      const VPM& point_map,
      const TM& tmesh, 
      typename boost::graph_traits<TM>::face_descriptor f, 
      const GeomTraits& gt) 
  {      
    // initialize all necessary kernel functions
    auto unit_normal = gt.construct_unit_normal_3_object();

    // reference and descriptor types
    typedef typename boost::property_traits<VPM>::reference Point_reference;
    typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;

    const halfedge_descriptor h = halfedge(f, tmesh);

    // get the three points of the face and calculate their unit normal
    const Point_reference p = get(point_map, source(h, tmesh));
    const Point_reference q = get(point_map, target(h, tmesh));
    const Point_reference r = get(point_map, target(next(h, tmesh), tmesh));

    return unit_normal(p, q, r);
  }
} //namespace common
} //namespace internal
} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_INTERNAL_COMMON_H
