// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Default.h>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

/*
This policy is not useful on its own; it is designed to be combined with another policy using a small weight.
Therefore, it is kept internal.
*/
template <typename TriangleMesh, typename GeomTraits>
class Line_quadric_calculator
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Row_4             Row_4;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

private:
  std::function<Vector_3(vertex_descriptor)> m_vertex_normal_map;

public:
  Line_quadric_calculator() = delete;

  template <typename VNM>
  Line_quadric_calculator(const VNM vnm)
      : m_vertex_normal_map([vnm](vertex_descriptor v) -> Vector_3{ return get(vnm, v); })
  { }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_vertex(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                                      const TriangleMesh& /*tmesh*/,
                                      const VertexPointMap point_map,
                                      const GeomTraits& gt) const
  {
    return construct_line_quadric_from_normal(m_vertex_normal_map(v), get(point_map, v), gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor /*he*/,
                                    const TriangleMesh& /*tmesh*/,
                                    const VertexPointMap /*point_map*/,
                                    const GeomTraits& /*gt*/) const
  {
    return Mat_4::Zero();
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor /*f*/,
                                    const TriangleMesh& /*tmesh*/,
                                    const VertexPointMap /*point_map*/,
                                    const GeomTraits& /*gt*/) const
  {
    return Mat_4::Zero();
  }

  /*
  // Since these policies are never used alone but always composed with another one, this code is not used
  Col_4 construct_optimal_point(const Mat_4& quadric,
                                const Col_4& p0,
                                const Col_4& p1) const
  {
    return construct_optimal_point_singular<GeomTraits>(quadric, p0, p1);
  }
  */
};

template<typename TriangleMesh,
         typename GeomTraits,
         typename VertexNormalMap = typename boost::property_map<TriangleMesh,
                                                CGAL::dynamic_vertex_property_t<typename GeomTraits::Vector_3> >::type >
class GarlandHeckbert_line_policies
  : public internal::GarlandHeckbert_cost_and_placement<
             internal::Line_quadric_calculator<TriangleMesh, GeomTraits>, TriangleMesh, GeomTraits>
{
public:
  typedef internal::Line_quadric_calculator<TriangleMesh, GeomTraits> Quadric_calculator;

private:
  typedef internal::GarlandHeckbert_cost_and_placement<
            Quadric_calculator, TriangleMesh, GeomTraits>                           Base;
  typedef GarlandHeckbert_line_policies<TriangleMesh, GeomTraits, VertexNormalMap>  Self;

public:
  typedef Self                                                                      Get_cost;
  typedef Self                                                                      Get_placement;

  typedef typename GeomTraits::FT                                                   FT;

public:
  template<typename VNM>
  GarlandHeckbert_line_policies(TriangleMesh& tmesh, const VNM vnm)
    : Base(tmesh, Quadric_calculator(vnm), 0)
  { }

public:
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  using Base::operator();
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H