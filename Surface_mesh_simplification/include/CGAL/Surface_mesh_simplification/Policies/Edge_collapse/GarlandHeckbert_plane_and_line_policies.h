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

namespace CGAL {
namespace Surface_mesh_simplification {

namespace internal{

// Helpers to perform computation in the constructor GH_plane_and_line_policies before calling the constructor of Base
template< typename TriangleMesh, typename NamedParameters>
struct GH_helper{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::FT FT;

  typedef dynamic_vertex_property_t<Vector_3>                                       Vertex_normal_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_normal_tag>::const_type Vertex_normal_dmap;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
                                                       NamedParameters,
                                                       Vertex_normal_dmap>::type    Vertex_normal_map;

  Vertex_normal_map vnm(const TriangleMesh& tmesh, const NamedParameters& np) const{
    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;
    constexpr bool must_compute_vertex_normals = is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value;

    Vertex_normal_map vertex_normals = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
                                                        get(Vertex_normal_tag(), tmesh));
    if constexpr(must_compute_vertex_normals)
      Polygon_mesh_processing::compute_vertex_normals(tmesh, vertex_normals, np);
    return vertex_normals;
  }

  double lw(const NamedParameters& np) const{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    return choose_parameter(get_parameter(np, internal_np::line_policies_weight), 0.01);
  }

  FT dm(const NamedParameters& np) const {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    // choose_parameter(get_parameter(np, internal_np::discontinuity_multiplier), 100);
    return FT(100);
  }
};

}

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

private:
  using TM=TriangleMesh;
  using GT=GeomTraits;


public:
  typedef typename Base::Quadric_calculator Quadric_calculator;

  typedef typename GeomTraits::FT                                              FT;

public:
  // GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
  //                                         const FT line_weight = FT(0.001),
  //                                         const FT dm = FT(100))
  //   : Base(tmesh, GH_plane_polices(tmesh, dm), GH_line_polices(tmesh, dm), FT(1.)/line_weight, dm)
  // { }

  template<typename VNM>
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
                                          const FT line_weight,
                                          const FT dm,
                                          const VNM vnm)
    : Base(tmesh, GH_plane_polices(tmesh, dm), GH_line_polices(tmesh, dm, vnm), FT(1.)/line_weight, dm)
  { }

  template<typename NP = parameters::Default_named_parameters>
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh, const NP& np = parameters::default_values()):
    Base(tmesh,
      GH_plane_polices(tmesh, internal::GH_helper<TM,NP>().dm(np)),
      GH_line_polices(tmesh, internal::GH_helper<TM,NP>().vnm(tmesh, np)),
      FT(1.)/internal::GH_helper<TM,NP>().lw(np),
      FT(1.),
      internal::GH_helper<TM,NP>().dm(np))
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