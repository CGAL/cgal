// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PMP_INTERNAL_CURVATURE_H
#define CGAL_PMP_INTERNAL_CURVATURE_H

#include <CGAL/license/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {

template <class TriangleMesh, class NamedParameters>
typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT
vertex_discrete_gaussian_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
                                    const TriangleMesh& tm,
                                    const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_const_property_map(vertex_point, tm));

  typename GeomTraits::Construct_vector_3 vector = gt.construct_vector_3_object();
  typename GeomTraits::Compute_scalar_product_3 scalar_product = gt.compute_scalar_product_3_object();
  typename GeomTraits::Compute_squared_length_3 squared_length = gt.compute_squared_length_3_object();
  typename GeomTraits::Construct_cross_product_vector_3 cross_product = gt.construct_cross_product_vector_3_object();

  const FT two_pi = 2 * CGAL_PI;

  FT ki = 0;
  for(halfedge_descriptor h : CGAL::halfedges_around_target(vd, tm))
  {
    if(is_border(h, tm))
      continue;

    const Vector_3 v0 = vector(get(vpm, vd), get(vpm, target(next(h, tm), tm))); // p1p2
    const Vector_3 v1 = vector(get(vpm, vd), get(vpm, source(h, tm))); // p1p0

    const FT dot = scalar_product(v0, v1);
    const Vector_3 cross = cross_product(v0, v1);
    const FT sqcn = squared_length(cross);
    if(dot == FT(0))
      ki += CGAL_PI/FT(2);
    else
    {
      if (sqcn == FT(0))
      {
        if (dot < 0)
          ki += CGAL_PI;
      }
      else{
        ki += std::atan2(CGAL::approximate_sqrt(sqcn), dot);
      }
    }
  }

  return two_pi - ki;
}

template <typename TriangleMesh, typename VertexCurvatureMap, typename NamedParameters>
void discrete_gaussian_curvature(const TriangleMesh& tm,
                                 VertexCurvatureMap vcm,
                                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;

  for(vertex_descriptor v : vertices(tm))
    put(vcm, v, vertex_discrete_gaussian_curvature(v, tm, np));
}


template <class TriangleMesh, class NamedParameters>
typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT
vertex_discrete_mean_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
                                    const TriangleMesh& tm,
                                    const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_const_property_map(vertex_point, tm));

  typename GeomTraits::Compute_squared_distance_3 squared_distance = gt.compute_squared_distance_3_object();
  typename GeomTraits::Compute_approximate_dihedral_angle_3 approximate_dihedral_angle= gt.compute_approximate_dihedral_angle_3_object();

  const FT two_pi = 2 * CGAL_PI;

  FT hi = 0;
  for(halfedge_descriptor h : CGAL::halfedges_around_target(vd, tm))
  {
    const Point_3& p = get(vpm, source(h, tm));
    const Point_3& q = get(vpm, target(h, tm));
    const Point_3& r = get(vpm, target(next(h, tm), tm));
    const Point_3& s = get(vpm, target(next(opposite(h, tm), tm), tm));
    const FT l = squared_distance(p,q);

    FT phi = CGAL_PI * approximate_dihedral_angle(p, q, r, s) / FT(180);

    if(phi < 0)
      phi += two_pi;
    if(phi > two_pi)
      phi = two_pi;

    hi += FT(0.5) * l * (CGAL_PI - phi);
  }

  return FT(0.5) * hi;
}

template <typename TriangleMesh, typename VertexCurvatureMap, typename NamedParameters>
void discrete_mean_curvature(const TriangleMesh& tm,
                                 VertexCurvatureMap vcm,
                                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;

  for(vertex_descriptor v : vertices(tm))
    put(vcm, v, vertex_discrete_mean_curvature(v, tm, np));
}

/// convenience overloads

template <class TriangleMesh>
auto
vertex_discrete_gaussian_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
                                    const TriangleMesh& tm)
{
  return vertex_discrete_gaussian_curvature(vd, tm, parameters::all_default());
}

template <typename TriangleMesh, typename VertexCurvatureMap>
void discrete_gaussian_curvature(const TriangleMesh& tm,
                                 VertexCurvatureMap vcm)
{
  discrete_gaussian_curvature(tm, vcm, parameters::all_default());
}

template <class TriangleMesh>
auto
vertex_discrete_mean_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
                                    const TriangleMesh& tm)
{
  return vertex_discrete_mean_curvature(vd, tm, parameters::all_default());
}

template <typename TriangleMesh, typename VertexCurvatureMap>
void discrete_mean_curvature(const TriangleMesh& tm,
                                 VertexCurvatureMap vcm)
{
  discrete_mean_curvature(tm, vcm, parameters::all_default());
}

} } // CGAL::Polygon_mesh_processing::experimental

#endif //CGAL_PMP_INTERNAL_CURVATURE_H
