// Copyright (c) 2018-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Konstantinos Katrioplas
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_OBB_H
#define CGAL_OPTIMAL_BOUNDING_BOX_OBB_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/population.h>
#include <CGAL/Optimal_bounding_box/internal/evolution.h>
#include <CGAL/Optimal_bounding_box/Optimal_bounding_box_traits.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
#include <CGAL/Real_timer.h>
#endif

#include <array>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

// works on matrices only
template <typename PointRange, typename Traits>
void post_processing(const typename Traits::Matrix& R,
                     std::array<typename Traits::Point_3, 8>& obb,
                     const PointRange& points,
                     const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Traits::Point_3                                   Point;
  typedef typename Traits::Matrix                                    Matrix;

  CGAL_assertion(R.number_of_rows() == 3 && R.number_of_columns() == 3);

  const Matrix Rt = traits.transpose(R);

  CGAL::Bbox_3 bbox;
  for(const Point& pt : points)
  {
    // @fixme should it be R here Rt at the other one... ?
    const FT x = pt.x(), y = pt.y(), z = pt.z();
    Point rotated_pt(x*R(0, 0) + y*R(0, 1) + z*R(0, 2),
                     x*R(1, 0) + y*R(1, 1) + z*R(1, 2),
                     x*R(2, 0) + y*R(2, 1) + z*R(2, 2));

    bbox += traits.construct_bbox_3_object()(rotated_pt);
  }

  // @todo could avoid building a cuboid
  typename Traits::Iso_cuboid_3 ic(bbox);

  // 3) apply inverse rotation to rotated AABB
  for(std::size_t i = 0; i<8; ++i)
  {
    const FT x = ic[i].x(), y = ic[i].y(), z = ic[i].z();
    obb[i] = Point(x*Rt(0, 0) + y*Rt(0, 1) + z*Rt(0, 2),
                   x*Rt(1, 0) + y*Rt(1, 1) + z*Rt(1, 2),
                   x*Rt(2, 0) + y*Rt(2, 1) + z*Rt(2, 2));
  }
}

template <typename PointRange, typename Traits>
void construct_optimal_bounding_box(std::array<typename Traits::Point_3, 8>& obb_points,
                                    const PointRange& points,
                                    CGAL::Random& rng,
                                    const Traits& traits)
{
  typedef typename Traits::Matrix                                    Matrix;

  const std::size_t max_generations = 100;

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  CGAL::Real_timer timer;
  timer.start();
#endif

  Evolution<PointRange, Traits> search_solution(points, rng, traits);
  search_solution.evolve(max_generations);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  std::cout << "evolve: " << timer.time() << std::endl;
  timer.reset();
#endif

  const Matrix& rotation = search_solution.get_best();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  std::cout << "get best: " << timer.time() << std::endl;
  timer.reset();
#endif

  post_processing(rotation, obb_points, points, traits);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  std::cout << "post-processing: " << timer.time() << std::endl;
#endif
}

template <typename PointRange, typename Traits>
void construct_optimal_bounding_box(std::array<typename Traits::Point_3, 8>& obb_points,
                                    const bool use_ch,
                                    const PointRange& points,
                                    CGAL::Random& rng,
                                    const Traits& traits)
{
  typedef typename Traits::Matrix                                    Matrix;
  typedef typename Traits::Point_3                                   Point;

  CGAL_static_assertion((std::is_same<typename boost::range_value<PointRange>::type, Point>::value));

  if(use_ch) // construct the convex hull to reduce the number of points
  {
    std::vector<Point> ch_points;
    extreme_points_3(points, std::back_inserter(ch_points));
    return construct_optimal_bounding_box(obb_points, ch_points, rng, traits);
  }
  else
  {
    return construct_optimal_bounding_box(obb_points, points, rng, traits);
  }
}

} // namespace Optimal_bounding_box
} // namespace internal

/// \ingroup PkgOptimalBoundingBoxFunctions
///
/// constructs a rectangular box that contains all the input points. This bounding box
/// is obtained via an optimization process aiming to get a close approximation of the
/// optimal bounding box, which is defined as the smallest (in terms of volume)
/// of all the rectangular boxes containing the input points.
///
/// \tparam PointRange a model of `Range` with value type `Point`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param points the input points
/// \param obb_points the eight points of the resulting englobing box.
///                   The order of points is the same as in the function `CGAL::make_hexahedron()`
/// \param np an optional sequence of \ref obb_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `pmesh`.
///     If this parameter is omitted, an internal property map for
///     `CGAL::vertex_point_t` must be available in `PolygonMesh`
///   \cgalParamEnd
///   \cgalParamBegin{geom_traits}
///     a geometric traits class instance, model of the concept `OptimalBoundingBoxTraits`. %Default is
///     `CGAL::Optimal_bounding_box::Optimal_bounding_box_traits<K>` where `K` is deduced
///     from the point type, which must then be compatible with `CGAL::Kernel_traits`.
///   \cgalParamEnd
///   \cgalParamBegin{use_convex_hull}
///     a Boolean value to indicate whether the algorithm should first extract the so-called extreme
///     points of the data range (i.e. construct the convex hull) to reduce the input data range
///     and accelerate the algorithm. The optimal value of this parameter will depend on the data
///     as it is a balance between two costs. %Default is `true`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \pre the value type of `PointRange` is `Point`
///
template <typename PointRange,
          typename Point,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void optimal_bounding_box(const PointRange& points,
                          std::array<Point, 8>& obb_points,
                          const CGAL_BGL_NP_CLASS& np)
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

#if defined(CGAL_EIGEN3_ENABLED)
  typedef typename CGAL::Kernel_traits<Point>::type                                 K;
  typedef Optimal_bounding_box::Optimal_bounding_box_traits<K>                      Default_traits;
#else
  typedef void                                                                      Default_traits;
#endif

  typedef typename internal_np::Lookup_named_param_def<internal_np::face_size_map_t,
                                                       CGAL_BGL_NP_CLASS,
                                                       Default_traits>::type        Geom_traits;

  CGAL_static_assertion_msg(!(std::is_same<Geom_traits, void>::value),
                            "You must provide a traits class or have Eigen enabled!");

  Geom_traits traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Geom_traits());

  const bool use_ch = choose_parameter(get_parameter(np, internal_np::use_convex_hull), false);
  const unsigned int seed = choose_parameter(get_parameter(np, internal_np::random_seed), 0); // undocumented

  CGAL::Random rng(seed);

  // @todo handle those cases instead
  CGAL_assertion(points.size() >= 3);
  if(points.size() <= 3)
  {
    std::cerr << "The optimal bounding box cannot YET be computed for a mesh with fewer than 4 vertices!\n";
    return;
  }

  return Optimal_bounding_box::internal::construct_optimal_bounding_box(obb_points, use_ch, points, rng, traits);
}

/// \cond SKIP_IN_MANUAL

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Convenience overloads for point ranges
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PointRange, typename Point>
void optimal_bounding_box(const PointRange& points,
                          std::array<Point, 8>& obb_points)
{
  return optimal_bounding_box(points, obb_points, CGAL::parameters::all_default());
}

/// \endcond

/// \ingroup PkgOptimalBoundingBoxFunctions
///
/// constructs a rectangular box that contains the input mesh. This bounding box
/// is obtained via an optimization process aiming to get a close approximation of the
/// optimal bounding box, which is defined as the smallest (in terms of volume)
/// of all the rectangular boxes containing the input mesh.
///
/// \tparam PolygonMesh a model of `FaceListGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param pmesh the input mesh
/// \param obb_mesh the resulting enclosing bounding box (an hexahedron)
/// \param np an optional sequence of \ref obb_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `pmesh`.
///     If this parameter is omitted, an internal property map for
///     `CGAL::vertex_point_t` must be available in `PolygonMesh`
///   \cgalParamEnd
///   \cgalParamBegin{geom_traits}
///     a geometric traits class instance, model of the concept `OptimalBoundingBoxTraits`. %Default is
///     `CGAL::Optimal_bounding_box::Optimal_bounding_box_traits<K>` where `K` is deduced
///     from the point type, which must be compatible with `CGAL::Kernel_traits`.
///   \cgalParamEnd
///   \cgalParamBegin{use_convex_hull}
///     a Boolean value to indicate whether the algorithm should first extract the so-called extreme
///     points of the data range (i.e. construct the convex hull) to reduce the input data range
///     and accelerate the algorithm. The optimal value of this parameter will depend on the data
///     as it is a balance between two costs. %Default is `true`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PolygonMesh,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void optimal_bounding_box(const PolygonMesh& pmesh,
                          PolygonMesh& obb_mesh,
                          const CGAL_BGL_NP_CLASS& np)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                  vertex_descriptor;

  typedef typename PMP::GetVertexPointMap<PolygonMesh, CGAL_BGL_NP_CLASS>::const_type   VPM;
  typedef typename boost::property_traits<VPM>::value_type                              Point;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  std::vector<Point> points;
  points.reserve(num_vertices(pmesh));

  for(vertex_descriptor v : vertices(pmesh))
    points.push_back(get(vpm, v));

  std::array<Point, 8> obb_points;
  optimal_bounding_box(points, obb_points, np);

  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_mesh);
}

/// \cond SKIP_IN_MANUAL

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Convenience overloads for polygon meshes
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PolygonMesh>
void optimal_bounding_box(const PolygonMesh& pmesh,
                          PolygonMesh& obb_mesh)
{
  return optimal_bounding_box(pmesh, obb_mesh, CGAL::parameters::all_default());
}

/// \endcond

} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_OBB_H
