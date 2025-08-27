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
// Author(s)     : Andreas Fabri,
//                 Mael Rouxel-Labbé

#ifndef CGAL_PMP_CURVATURE_H
#define CGAL_PMP_CURVATURE_H

#include <CGAL/license/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Weights/cotangent_weights.h>

#include <cmath>
#include <algorithm>

namespace CGAL {
namespace Polygon_mesh_processing {

/**
  * \ingroup PMP_measure_grp
  *
  * computes the sum of the angles around a vertex.
  *
  * The angle sum is given in degrees.
  *
  * @tparam PolygonMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param v the vertex whose sum of angles is computed
  * @param pmesh the polygon mesh to which `v` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must provide the nested functor `Compute_approximate_angle_3`,
 *                    model of `Kernel::ComputeApproximateAngle_3`.}
 *     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the sum of angles around `v`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `pmesh`.
  *
  * \warning This function involves trigonometry.
  */
template<typename PolygonMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type::FT
#endif
angle_sum(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
          const PolygonMesh& pmesh,
          const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Geom_traits = typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;
  using FT = typename Geom_traits::FT;

  typename GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  CGAL_precondition(is_valid_vertex_descriptor(v, pmesh));

  typename Geom_traits::Compute_approximate_angle_3 approx_angle = gt.compute_approximate_angle_3_object();

  FT angle_sum = 0;
  for(auto h : halfedges_around_source(v, pmesh))
  {
    if(is_border(h, pmesh))
      continue;

    angle_sum += approx_angle(get(vpm, target(h, pmesh)),
                              get(vpm, source(h, pmesh)),
                              get(vpm, source(prev(h,pmesh), pmesh)));
  }

  return angle_sum;
}

// Discrete Gaussian Curvature

/**
  * \ingroup PMP_measure_grp
  *
  * computes the discrete Gaussian curvature at a vertex.
  *
  * We refer to Meyer et al. \cgalCite{cgal:mdsb-ddgot-02} for the definition of <i>discrete Gaussian curvature</i>.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param v the vertex whose discrete Gaussian curvature is being computed
  * @param tmesh the triangle mesh to which `v` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must be a model of `Kernel`.}
 *     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the discrete Gaussian curvature at `v`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * \warning This function involves trigonometry.
  * \warning The current formulation is not well defined for border vertices.
  *
  * \pre `tmesh` is a triangle mesh
  */
template <typename TriangleMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
discrete_Gaussian_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                            const TriangleMesh& tmesh,
                            const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GeomTraits = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type;
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  using VertexPointMap = typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_const_property_map(vertex_point, tmesh));

  typename GeomTraits::Construct_vector_3 vector =
    gt.construct_vector_3_object();
  typename GeomTraits::Construct_cross_product_vector_3 cross_product =
    gt.construct_cross_product_vector_3_object();
  typename GeomTraits::Compute_scalar_product_3 scalar_product =
    gt.compute_scalar_product_3_object();
  typename GeomTraits::Compute_squared_length_3 squared_length =
    gt.compute_squared_length_3_object();

  FT angle_sum = 0;

  for(halfedge_descriptor h : CGAL::halfedges_around_target(v, tmesh))
  {
    if(is_border(h, tmesh))
      continue;

    const Vector_3 v0 = vector(get(vpm, v), get(vpm, target(next(h, tmesh), tmesh))); // p1p2
    const Vector_3 v1 = vector(get(vpm, v), get(vpm, source(h, tmesh))); // p1p0

    const FT dot = scalar_product(v0, v1);
    const Vector_3 cross = cross_product(v0, v1);
    const FT sqcn = squared_length(cross);
    if(is_zero(dot))
    {
      angle_sum += CGAL_PI/FT(2);
    }
    else
    {
      if(is_zero(sqcn)) // collinear
      {
        if(dot < 0)
          angle_sum += CGAL_PI;
        // else
        //   angle_sum += 0;
      }
      else
      {
        angle_sum += std::atan2(CGAL::approximate_sqrt(sqcn), dot);
      }
    }
  }

  Weights::Secure_cotangent_weight_with_voronoi_area<TriangleMesh, VertexPointMap, GeomTraits> wc(tmesh, vpm, gt);

  const FT gaussian_curvature = (2 * CGAL_PI - angle_sum) / wc.voronoi(v);

  return gaussian_curvature;
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the discrete Gaussian curvatures at the vertices of a mesh.
  *
  * We refer to Meyer et al. \cgalCite{cgal:mdsb-ddgot-02} for the definition of <i>discrete Gaussian curvature</i>.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam VertexCurvatureMap must be a model of `WritablePropertyMap` with key type
  *                            `boost::graph_traits<TriangleMesh>::%vertex_descriptor` and value type `FT`,
  *                            which is either `geom_traits::FT` if this named parameter is provided,
  *                            or `kernel::FT` with the kernel deduced from from the point property map of `tmesh`.
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tmesh the triangle mesh to which `v` belongs
  * @param vcm the property map that contains the computed discrete curvatures
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must be a model of `Kernel`.}
 *     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * \warning This function involves trigonometry.
  * \warning The current formulation is not well defined for border vertices.
  *
  * \pre `tmesh` is a triangle mesh
  */
template <typename TriangleMesh,
          typename VertexCurvatureMap,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void discrete_Gaussian_curvatures(const TriangleMesh& tmesh,
                                  VertexCurvatureMap vcm,
                                  const CGAL_NP_CLASS& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  for(vertex_descriptor v : vertices(tmesh))
  {
    put(vcm, v, discrete_Gaussian_curvature(v, tmesh, np));
    // std::cout << "curvature: " << get(vcm, v) << std::endl;
  }
}

// Discrete Mean Curvature

/**
  * \ingroup PMP_measure_grp
  *
  * computes the discrete mean curvature at a vertex.
  *
  * We refer to Meyer et al. \cgalCite{cgal:mdsb-ddgot-02} for the definition of <i>discrete mean curvature</i>.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param v the vertex whose discrete mean curvature is being computed
  * @param tmesh the triangle mesh to which `v` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must be a model of `Kernel`.}
 *     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the discrete mean curvature at `v`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * \warning The current formulation is not well defined for border vertices.
  *
  * \pre `tmesh` is a triangle mesh
  */
template <typename TriangleMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
discrete_mean_curvature(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                        const TriangleMesh& tmesh,
                        const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GeomTraits = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type;
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  using VertexPointMap = typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type;
  using Point_ref = typename boost::property_traits<VertexPointMap>::reference;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_const_property_map(vertex_point, tmesh));

#if 0
  typename GeomTraits::Compute_squared_distance_3 squared_distance =
    gt.compute_squared_distance_3_object();
  typename GeomTraits::Compute_approximate_dihedral_angle_3 dihedral_angle =
    gt.compute_approximate_dihedral_angle_3_object();

  const FT two_pi = 2 * CGAL_PI;

  FT hi = 0;
  for(halfedge_descriptor h : CGAL::halfedges_around_target(v, tmesh))
  {
    const Point_3& p = get(vpm, source(h, tmesh));
    const Point_3& q = get(vpm, target(h, tmesh));
    const Point_3& r = get(vpm, target(next(h, tmesh), tmesh));
    const Point_3& s = get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh));
    const FT l = squared_distance(p,q);

    FT phi = CGAL_PI * dihedral_angle(p, q, r, s) / FT(180);

    if(phi < 0)
      phi += two_pi;
    if(phi > two_pi)
      phi = two_pi;

    hi += FT(0.5) * l * (CGAL_PI - phi);
  }

  return FT(0.5) * hi;
#else
  typename GeomTraits::Construct_vector_3 vector =
    gt.construct_vector_3_object();
  typename GeomTraits::Construct_sum_of_vectors_3 vector_sum =
    gt.construct_sum_of_vectors_3_object();
  typename GeomTraits::Construct_scaled_vector_3 scaled_vector =
    gt.construct_scaled_vector_3_object();
  typename GeomTraits::Compute_squared_length_3 squared_length =
    gt.compute_squared_length_3_object();

  Weights::Secure_cotangent_weight_with_voronoi_area<TriangleMesh, VertexPointMap, GeomTraits> wc(tmesh, vpm, gt);

  Vector_3 kh = vector(CGAL::NULL_VECTOR);
  for(halfedge_descriptor h : CGAL::halfedges_around_target(v, tmesh))
  {
    const vertex_descriptor v1 = source(h, tmesh);

    const Point_ref p0 = get(vpm, v);
    const Point_ref p1 = get(vpm, v1);

    FT local_c = 0;
    if(!is_border(h, tmesh))
    {
      const vertex_descriptor v2 = target(next(h, tmesh), tmesh);
      const Point_ref p2 = get(vpm, v2);
      local_c += Weights::cotangent_3_clamped(p0, p2, p1, gt);
    }

    if(!is_border(opposite(h, tmesh), tmesh))
    {
      const vertex_descriptor v3 = target(next(opposite(h, tmesh), tmesh), tmesh);
      const Point_ref p3 = get(vpm, v3);
      local_c += Weights::cotangent_3_clamped(p1, p3, p0, gt);
    }

    kh = vector_sum(kh, scaled_vector(vector(p0, p1), local_c));
  }

  const FT khn = CGAL::approximate_sqrt(squared_length(kh));
  const FT va = wc.voronoi(v);
  CGAL_assertion(!is_zero(va));

  const FT mean_curvature = khn / (FT(4) * va);
  return mean_curvature;
#endif
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the discrete mean curvatures at the vertices of a mesh.
  *
  * We refer to Meyer et al. \cgalCite{cgal:mdsb-ddgot-02} for the definition of <i>discrete mean curvature</i>.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam VertexCurvatureMap must be a model of `WritablePropertyMap` with key type
  *                            `boost::graph_traits<TriangleMesh>::%vertex_descriptor` and value type `FT`,
  *                            which is either `geom_traits::FT` if this named parameter is provided,
  *                            or `kernel::FT` with the kernel deduced from from the point property map of `tmesh`.
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tmesh the triangle mesh to which `v` belongs
  * @param vcm the property map that contains the computed discrete curvatures
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must be a model of `Kernel`.}
 *     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * \warning The current formulation is not well defined for border vertices.
  *
  * \pre `tmesh` is a triangle mesh
  */
template <typename TriangleMesh,
          typename VertexCurvatureMap,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void discrete_mean_curvatures(const TriangleMesh& tmesh,
                              VertexCurvatureMap vcm,
                              const CGAL_NP_CLASS& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  for(vertex_descriptor v : vertices(tmesh))
    put(vcm, v, discrete_mean_curvature(v, tmesh, np));
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //CGAL_PMP_CURVATURE_H
