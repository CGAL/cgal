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

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template <typename TriangleMesh, typename GeomTraits>
class Plane_quadric_calculator
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

public:
  Plane_quadric_calculator() { }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return construct_classic_plane_quadric_from_edge(he, tmesh, point_map, gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return construct_classic_plane_quadric_from_face(f, tmesh, point_map, gt);
  }

  // @fixme unused?
  Mat_4 construct_quadric_from_normal(const typename GeomTraits::Vector_3& normal,
                                      const typename GeomTraits::Point_3& point,
                                      const GeomTraits& gt) const
  {
    typedef typename GeomTraits::FT                                            FT;

    auto dot_product = gt.compute_scalar_product_3_object();
    auto construct_vector = gt.construct_vector_3_object();

    // negative dot product between the normal and the position vector
    const FT d = - dot_product(normal, construct_vector(ORIGIN, point));

    // row vector given by d appended to the normal
    const Col_4 row { normal.x(), normal.y(), normal.z(), d };

    // outer product
    return row.transpose() * row;
  }

  Col_4 construct_optimal_point(const Mat_4& quadric,
                                const Col_4& p0,
                                const Col_4& p1) const
  {
    return construct_optimal_point_singular<GeomTraits>(quadric, p0, p1);
  }
};

} // namespace internal

template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_policies
  : public internal::GarlandHeckbert_cost_and_placement<
             internal::Plane_quadric_calculator<TriangleMesh, GeomTraits>, TriangleMesh, GeomTraits>
{
public:
  typedef internal::Plane_quadric_calculator<TriangleMesh, GeomTraits>         Quadric_calculator;

private:
  typedef internal::GarlandHeckbert_cost_and_placement<
            Quadric_calculator, TriangleMesh, GeomTraits>                      Base;
  typedef GarlandHeckbert_plane_policies<TriangleMesh, GeomTraits>             Self;

public:
  typedef Self                                                                 Get_cost;
  typedef Self                                                                 Get_placement;

  typedef typename GeomTraits::FT                                              FT;

public:
  GarlandHeckbert_plane_policies(TriangleMesh& tmesh,
                                 const FT dm = FT(100))
    : Base(tmesh, Quadric_calculator(), dm)
  { }

public:
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  using Base::operator();
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_POLICIES_H
