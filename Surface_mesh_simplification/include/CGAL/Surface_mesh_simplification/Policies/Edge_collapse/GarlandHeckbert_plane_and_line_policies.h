// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_AND_LINE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_AND_LINE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_composed_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_line_policies.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_plane_policies.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<typename TriangleMesh,
         typename GeomTraits,
         typename VertexNormalMap = typename boost::property_map<TriangleMesh,
                                                CGAL::dynamic_vertex_property_t<typename GeomTraits::Vector_3> >::const_type >
class GarlandHeckbert_plane_and_line_policies
  : public internal::GarlandHeckbert_composed_policies<TriangleMesh, GeomTraits,
                                             GarlandHeckbert_plane_policies<TriangleMesh, GeomTraits>,
                                             internal::GarlandHeckbert_line_policies<TriangleMesh, GeomTraits, VertexNormalMap>,
                                             true>
{
  typedef GarlandHeckbert_plane_policies<TriangleMesh, GeomTraits>                            GH_plane_polices;
  typedef internal::GarlandHeckbert_line_policies<TriangleMesh, GeomTraits, VertexNormalMap>  GH_line_polices;
  typedef internal::GarlandHeckbert_composed_policies<TriangleMesh, GeomTraits,
                                                      GH_plane_polices,
                                                      GH_line_polices,
                                                      true>                                   Base;


public:
  typedef typename Base::Quadric_calculator Quadric_calculator;

  typedef typename GeomTraits::FT                                              FT;

public:
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
                                          const FT line_weight = FT(0.01),
                                          const FT dm = FT(100))
    : Base(tmesh, GH_plane_polices(tmesh, dm), GH_line_polices(tmesh, dm), FT(1.)/line_weight, dm)
  { }

  template<typename VNM>
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
                                          const FT line_weight,
                                          const FT dm,
                                          const VNM vnm)
    : Base(tmesh, GH_plane_polices(tmesh, dm), GH_line_polices(tmesh, dm, vnm), FT(1.)/line_weight, dm)
  { }

public:
  using Base::operator();
  using Base::get_cost;
  using Base::get_placement;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_AND_LINE_POLICIES_H