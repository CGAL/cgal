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
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_ORIENTED_BOUNDING_BOX_H
#define CGAL_OPTIMAL_BOUNDING_BOX_ORIENTED_BOUNDING_BOX_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/population.h>
#include <CGAL/Optimal_bounding_box/internal/evolution.h>
#include <CGAL/Optimal_bounding_box/Oriented_bounding_box_traits.h>

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

template <typename PointRange, typename Traits>
void construct_oriented_bounding_box(std::array<typename Traits::Point_3, 8>& obb_points,
                                     const typename Traits::Aff_transformation_3& transformation,
                                     const typename Traits::Aff_transformation_3& inverse_transformation,
                                     const PointRange& points,
                                     const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Traits::Point_3                                   Point;

  // Construct the bbox of the transformed point set
  CGAL::Bbox_3 bbox;
  for(const Point& pt : points)
  {
    const Point rotated_pt = transformation.transform(pt);
    bbox += traits.construct_bbox_3_object()(rotated_pt);
  }

  obb_points[0] = Point(bbox.xmin(), bbox.ymin(), bbox.zmin());
  obb_points[1] = Point(bbox.xmax(), bbox.ymin(), bbox.zmin());
  obb_points[2] = Point(bbox.xmax(), bbox.ymax(), bbox.zmin());
  obb_points[3] = Point(bbox.xmin(), bbox.ymax(), bbox.zmin());

  obb_points[4] = Point(bbox.xmin(), bbox.ymax(), bbox.zmax()); // see order in make_hexahedron()...
  obb_points[5] = Point(bbox.xmin(), bbox.ymin(), bbox.zmax());
  obb_points[6] = Point(bbox.xmax(), bbox.ymin(), bbox.zmax());
  obb_points[7] = Point(bbox.xmax(), bbox.ymax(), bbox.zmax());

  // Apply the inverse rotation to the rotated axis aligned bounding box
  for(std::size_t i=0; i<8; ++i)
    obb_points[i] = inverse_transformation.transform(obb_points[i]);
}

template <typename PointRange, typename Traits>
void compute_best_transformation(typename Traits::Aff_transformation_3& transformation,
                                 typename Traits::Aff_transformation_3& inverse_transformation,
                                 const PointRange& points,
                                 CGAL::Random& rng,
                                 const Traits& traits)
{
  typedef typename Traits::Matrix                                    Matrix;
  typedef typename Traits::Aff_transformation_3                      Aff_transformation_3;

  const std::size_t max_generations = 50; // @todo hidden NP

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

  const Matrix& rot = search_solution.get_best();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  std::cout << "get best: " << timer.time() << std::endl;
#endif

  transformation = Aff_transformation_3(rot(0, 0), rot(0, 1), rot(0, 2),
                                        rot(1, 0), rot(1, 1), rot(1, 2),
                                        rot(2, 0), rot(2, 1), rot(2, 2));

  // inverse transformation is simply the transposed since the matrix is unitary
  inverse_transformation = Aff_transformation_3(rot(0, 0), rot(1, 0), rot(2, 0),
                                                rot(0, 1), rot(1, 1), rot(2, 1),
                                                rot(0, 2), rot(1, 2), rot(2, 2));
}

// Following two functions are overload to dispatch depending on return type
template <typename PointRange, typename Traits>
void construct_oriented_bounding_box(typename Traits::Aff_transformation_3& transformation,
                                     const PointRange& points,
                                     CGAL::Random& rng,
                                     const Traits& traits)
{
  typename Traits::Aff_transformation_3 inverse_transformation;
  compute_best_transformation(transformation, inverse_transformation, points, rng, traits);
}

template <typename PointRange, typename Traits>
void construct_oriented_bounding_box(std::array<typename Traits::Point_3, 8>& obb_points,
                                     const PointRange& points,
                                     CGAL::Random& rng,
                                     const Traits& traits)
{
  typename Traits::Aff_transformation_3 transformation, inverse_transformation;
  compute_best_transformation(transformation, inverse_transformation, points, rng, traits);

  construct_oriented_bounding_box(obb_points, transformation, inverse_transformation, points, traits);
}

// Entry point, decide whether to compute the CH_3 or not
template <typename Output, typename PointRange, typename Traits>
void construct_oriented_bounding_box(Output& output,
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
    std::cout << "points on CH: " << ch_points.size() << std::endl;
    return construct_oriented_bounding_box(output, ch_points, rng, traits);
  }
  else
  {
    return construct_oriented_bounding_box(output, points, rng, traits);
  }
}

} // namespace Optimal_bounding_box
} // namespace internal

/// \addtogroup PkgOptimalBoundingBox_Oriented_bounding_box
///
/// The function `oriented_bounding_box` computes an approximation of the <i>optimal bounding box</i>,
/// which is defined as the rectangular box with smallest volume of all the rectangular boxes containing
/// the input points.
///
/// Internally, the algorithm uses an optimization process to compute a transformation (rotation)
/// \f$ {\mathcal R}_b\f$ such that the axis-aligned box of the rotated input point set
/// has a volume that is as small as possible given a fixed maximal number of optimization iterations.
///
/// \cgalHeading{Input}
///
/// The input can be either a range of points, or a polygon mesh.
///
/// \cgalHeading{Output}
///
/// The result of the algorithm can be retrieved as either:
/// - the best affine transformation (\f$ {\mathcal R}_b\f$) that the algorithm has found;
/// - an array of eight points, representing the best oriented bounding box (\f$ {\mathcal B}_b\f$)
///   that the algorithm has constructed, which is related to (\f$ {\mathcal R}_b\f$) as it is
///   the inverse transformation of the axis-aligned bounding box of the transformed point set.
///   The order of the points in the array is the same as in the function
///   \link PkgBGLHelperFct `CGAL::make_hexahedron()` \endlink,
///   which is a useful function to construct a mesh from these points.
///
/// Note that when returning an array of points, these points are constructed from the axis-aligned
/// bounding box and some precision loss should therefore be expected if a kernel not providing
/// exact constructions is used.
///
/// The algorithm is based on a paper by Chang, Gorissen, and Melchior \cgalCite{cgal:cgm-fobbo-11}.

/// \ingroup PkgOptimalBoundingBox_Oriented_bounding_box
///
/// See above.
///
/// \tparam PointRange a model of `Range`
/// \tparam Output either `std::array<Point, 8>` with `Point` being equivalent to the traits' `Point_3` type,
///                or the traits' `Aff_transformation_3` type
/// \tparam NamedParameters a sequence of \ref obb_namedparameters "Named Parameters"
///
/// \param points the input points
/// \param out the resulting array of points or affine transformation
/// \param np an optional sequence of \ref obb_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `pmesh`.
///     If this parameter is omitted, an internal property map for
///     `CGAL::vertex_point_t` must be available in `PolygonMesh`
///   \cgalParamEnd
///   \cgalParamBegin{geom_traits}
///     a geometric traits class instance, model of the concept `OrientedBoundingBoxTraits`.
///     %Default is `CGAL::Oriented_bounding_box_traits<K>` where `K` is deduced from the point type.
///   \cgalParamEnd
///   \cgalParamBegin{use_convex_hull}
///     a Boolean value to indicate whether the algorithm should first extract the so-called extreme
///     points of the data range (i.e. construct the convex hull) to reduce the input data range
///     and accelerate the algorithm. %Default is `true`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PointRange,
          typename Output,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void oriented_bounding_box(const PointRange& points,
                           Output& out,
                           const CGAL_BGL_NP_CLASS& np,
#ifndef DOXYGEN_RUNNING
                           typename boost::enable_if<
                             typename boost::has_range_iterator<PointRange>
                                    >::type* = 0
#endif
                           )
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

#if defined(CGAL_EIGEN3_ENABLED)
  typedef typename boost::range_value<PointRange>::type                                 Point;
  typedef typename CGAL::Kernel_traits<Point>::type                                     K;
  typedef Oriented_bounding_box_traits<K>                                               Default_traits;
#else
  typedef void                                                                          Default_traits;
#endif

  typedef typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                       CGAL_BGL_NP_CLASS,
                                                       Default_traits>::type        Geom_traits;

  CGAL_static_assertion_msg(!(std::is_same<Geom_traits, void>::value),
                            "You must provide a traits class or have Eigen enabled!");

  Geom_traits traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Geom_traits());

  const bool use_ch = choose_parameter(get_parameter(np, internal_np::use_convex_hull), true);
  const unsigned int seed = choose_parameter(get_parameter(np, internal_np::random_seed), 0); // undocumented

  CGAL::Random rng(seed);

  // @todo handle those cases instead
  CGAL_assertion(points.size() >= 3);
  if(points.size() <= 3)
  {
    std::cerr << "The oriented bounding box cannot YET be computed for a mesh with fewer than 4 vertices!\n";
    return;
  }

  return Optimal_bounding_box::internal::construct_oriented_bounding_box(out, use_ch, points, rng, traits);
}

/// \ingroup PkgOptimalBoundingBox_Oriented_bounding_box
///
/// Extracts the vertices of the mesh as a point range and calls the other overload.
///
/// \tparam PolygonMesh a model of `VertexListGraph`
/// \tparam Output either `std::array<Point, 8>` with `Point` being equivalent to the traits' `Point_3` type,
///                or the traits' `Aff_transformation_3` type
/// \tparam NamedParameters a sequence of \ref obb_namedparameters "Named Parameters"
///
/// \param pmesh the input mesh
/// \param out the resulting array of points or affine transformation
/// \param np an optional sequence of \ref obb_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `pmesh`.
///     If this parameter is omitted, an internal property map for
///     `CGAL::vertex_point_t` must be available in `PolygonMesh`
///   \cgalParamEnd
///   \cgalParamBegin{geom_traits}
///     a geometric traits class instance, model of the concept `OrientedBoundingBoxTraits`.
///     %Default is `CGAL::Oriented_bounding_box_traits<K>` where `K` is deduced from the point type.
///   \cgalParamEnd
///   \cgalParamBegin{use_convex_hull}
///     a Boolean value to indicate whether the algorithm should first extract the so-called extreme
///     points of the data range (i.e. construct the convex hull) to reduce the input data range
///     and accelerate the algorithm. %Default is `true`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PolygonMesh,
          typename Output,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void oriented_bounding_box(const PolygonMesh& pmesh,
                           Output& out,
                           const CGAL_BGL_NP_CLASS& np,
#ifndef DOXYGEN_RUNNING
                           typename boost::disable_if<
                             typename boost::has_range_iterator<PolygonMesh>
                                    >::type* = 0
#endif
                           )
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

  oriented_bounding_box(points, out, np);
}

/// \cond SKIP_IN_MANUAL

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Convenience overloads
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename InputData /*range or mesh*/, typename OutputType /*array or transformation*/>
void oriented_bounding_box(const InputData& data,
                           OutputType& out)
{
  return oriented_bounding_box(data, out, CGAL::parameters::all_default());
}

/// \endcond

} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_ORIENTED_BOUNDING_BOX_H
