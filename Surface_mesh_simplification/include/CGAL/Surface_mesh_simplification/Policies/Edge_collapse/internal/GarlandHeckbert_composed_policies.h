// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_COMPOSED_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_COMPOSED_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

/*
That policy was created to merge the line policy with the existing ones.
Since it is a technical feature, it is internal to the library.
The composition of some quadrics is always invertible, but this is not necessarily true for all composed quadrics.
`invertible` boolean parameter reflects this distinction.

*/

template <typename TriangleMesh, typename GeomTraits, typename Quadric_calculator_1, typename Quadric_calculator_2, bool invertible=false>
class Composed_quadric_calculator
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Row_4             Row_4;

  Quadric_calculator_1 quadric_calculator_1;
  Quadric_calculator_2 quadric_calculator_2;
  double weight_1;
  double weight_2;

public:
  Composed_quadric_calculator(const Quadric_calculator_1& qc1, const Quadric_calculator_2& qc2, double w1 = 1., double w2 = 1.)
      : quadric_calculator_1(qc1)
      , quadric_calculator_2(qc2)
      , weight_1(w1)
      , weight_2(w2) {}
  Composed_quadric_calculator(double w1 = 1., double w2 = 1.)
      : weight_1(w1)
      , weight_2(w2) {}

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_vertex(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                                      const TriangleMesh& tmesh,
                                      const VertexPointMap point_map,
                                      const GeomTraits& gt) const
  {
    return weight_1 * quadric_calculator_1.construct_quadric_from_vertex(v, tmesh, point_map, gt) +
           weight_2 * quadric_calculator_2.construct_quadric_from_vertex(v, tmesh, point_map, gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return weight_1 * quadric_calculator_1.construct_quadric_from_edge(he, tmesh, point_map, gt) +
           weight_2 * quadric_calculator_2.construct_quadric_from_edge(he, tmesh, point_map, gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return weight_1 * quadric_calculator_1.construct_quadric_from_face(f, tmesh, point_map, gt) +
           weight_2 * quadric_calculator_2.construct_quadric_from_face(f, tmesh, point_map, gt);
  }


  Col_4 construct_optimal_point(const Mat_4& quadric,
                                const Col_4& p0,
                                const Col_4& p1) const
  {
    if constexpr(invertible)
      return construct_optimal_point_invertible<GeomTraits>(quadric);
    else
      return construct_optimal_point_singular<GeomTraits>(quadric, p0, p1);
  }
};

template<typename TriangleMesh, typename GeomTraits, typename GH_policies_1, typename GH_policies_2, bool invertible=false>
class GarlandHeckbert_composed_policies
  : public internal::GarlandHeckbert_cost_and_placement<
             internal::Composed_quadric_calculator<TriangleMesh, GeomTraits,
                                                   typename GH_policies_1::Quadric_calculator,
                                                   typename GH_policies_2::Quadric_calculator,
                                                   invertible>,
             TriangleMesh, GeomTraits>
{
public:
  typedef internal::Composed_quadric_calculator<TriangleMesh, GeomTraits,
                                               typename GH_policies_1::Quadric_calculator,
                                               typename GH_policies_2::Quadric_calculator,
                                               invertible>
          Quadric_calculator;

private:
  typedef internal::GarlandHeckbert_cost_and_placement<
            Quadric_calculator, TriangleMesh, GeomTraits>                      Base;
  typedef GarlandHeckbert_composed_policies<TriangleMesh, GeomTraits, GH_policies_1, GH_policies_2, invertible> Self;

public:
  typedef Self                                                                 Get_cost;
  typedef Self                                                                 Get_placement;

  typedef typename GeomTraits::FT                                              FT;

public:
  GarlandHeckbert_composed_policies(TriangleMesh& tmesh,
                                    double w1=1., double w2=1.,const FT dm = FT(100))
    : Base(tmesh, Quadric_calculator(w1, w2), dm)
  { }

  GarlandHeckbert_composed_policies(TriangleMesh& tmesh,
                                    GH_policies_1 ghp1,
                                    GH_policies_2 ghp2,
                                    double w1=1., double w2=1.,const FT dm = FT(100))
    : Base(tmesh, Quadric_calculator(ghp1.quadric_calculator(), ghp2.quadric_calculator(), w1, w2), dm)
  { }

public:
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  using Base::operator();
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_COMPOSED_POLICIES_H