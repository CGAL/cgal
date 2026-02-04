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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <boost/property_map/function_property_map.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

namespace internal {

// Helpers to perform computation in the constructor GH_plane_and_line_policies before calling the constructor of Base
template< typename TriangleMesh = void, typename NamedParameters = void>
struct GH_helper{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::FT FT;

  typedef dynamic_vertex_property_t<Vector_3>                                       Vertex_normal_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_normal_tag>::type       Vertex_normal_dmap;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
                                                       NamedParameters,
                                                       Vertex_normal_dmap>::type    Vertex_normal_map;

  NamedParameters np;
  GH_helper(const NamedParameters &np_):np(np_){ }

  Vertex_normal_map vnm(TriangleMesh& tmesh) const{
    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;
    constexpr bool must_compute_vertex_normals = is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value;

    Vertex_normal_map vertex_normals = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
                                                        get(Vertex_normal_tag(), tmesh));
    if constexpr(must_compute_vertex_normals){
      Polygon_mesh_processing::compute_vertex_normals(tmesh, vertex_normals, np);
    }
    return vertex_normals;
  }

  double lw() const{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    return choose_parameter(get_parameter(np, internal_np::line_policies_weight), 0.01);
  }

  double dm() const{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    return choose_parameter(get_parameter(np, internal_np::discontinuity_multiplier), 100);
  }
};

} // namespace internal

template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_and_line_policies
  : public internal::GarlandHeckbert_composed_policies<TriangleMesh, GeomTraits,
                                             GarlandHeckbert_plane_policies<TriangleMesh, GeomTraits>,
                                             internal::GarlandHeckbert_line_policies<TriangleMesh, GeomTraits>,
                                             true>
{
  typedef GarlandHeckbert_plane_policies<TriangleMesh, GeomTraits>                            GH_plane_polices;
  typedef internal::GarlandHeckbert_line_policies<TriangleMesh, GeomTraits>                   GH_line_polices;
  typedef internal::GarlandHeckbert_composed_policies<TriangleMesh, GeomTraits,
                                                      GH_plane_polices,
                                                      GH_line_polices,
                                                      true>                                   Base;

public:
  typedef typename Base::Quadric_calculator Quadric_calculator;

  typedef typename GeomTraits::FT                                              FT;

public:
  template<typename NP = parameters::Default_named_parameters>
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh, const NP& np = parameters::default_values()):
    Base(tmesh,
      GH_plane_polices(tmesh),
      GH_line_polices(tmesh, internal::GH_helper<TriangleMesh,NP>(np).vnm(tmesh)),
      FT(1.)/internal::GH_helper<TriangleMesh,NP>(np).lw(),
      FT(1.),
      internal::GH_helper<TriangleMesh,NP>(np).dm())
  { }

public:
  using Base::operator();
  using Base::get_cost;
  using Base::get_placement;
};

template<typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
GarlandHeckbert_plane_and_line_policies<TriangleMesh, typename GetGeomTraits<TriangleMesh,NamedParameters>::type>
  make_GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
                                               const NamedParameters& np = parameters::default_values()){
  typedef typename GetGeomTraits<TriangleMesh,NamedParameters>::type GT;
  return GarlandHeckbert_plane_and_line_policies<TriangleMesh, GT>(tmesh, np);
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_AND_LINE_POLICIES_H