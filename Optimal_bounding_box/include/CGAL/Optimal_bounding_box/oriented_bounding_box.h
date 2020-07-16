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

#include <CGAL/Optimal_bounding_box/internal/evolution.h>
#include <CGAL/Optimal_bounding_box/internal/population.h>
#include <CGAL/Optimal_bounding_box/Oriented_bounding_box_traits_3.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
#include <CGAL/Real_timer.h>
#endif

#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/value_type.hpp>
#include <boost/utility/enable_if.hpp>

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
void construct_oriented_bounding_box(const PointRange& points,
                                     const typename Traits::Aff_transformation_3& transformation,
                                     const typename Traits::Aff_transformation_3& inverse_transformation,
                                     std::array<typename Traits::Point_3, 8>& obb_points,
                                     const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Traits::Point_3                                   Point;

  CGAL_precondition(!points.empty());

  // Construct the bbox of the transformed point set
  typename PointRange::const_iterator pit = std::begin(points);
  const Point& first_pt = *pit++;
  const Point first_rot_pt = transformation.transform(first_pt);
  FT xmin = first_rot_pt.x(), xmax = first_rot_pt.x();
  FT ymin = first_rot_pt.y(), ymax = first_rot_pt.y();
  FT zmin = first_rot_pt.z(), zmax = first_rot_pt.z();

  for(typename PointRange::const_iterator end=std::end(points); pit!=end; ++pit)
  {
    const Point rot_pt = transformation.transform(*pit);

    xmin = (std::min)(rot_pt.x(), xmin);
    ymin = (std::min)(rot_pt.y(), ymin);
    zmin = (std::min)(rot_pt.z(), zmin);
    xmax = (std::max)(rot_pt.x(), xmax);
    ymax = (std::max)(rot_pt.y(), ymax);
    zmax = (std::max)(rot_pt.z(), zmax);
  }

  typename Traits::Construct_point_3 cp = traits.construct_point_3_object();

  obb_points[0] = cp(xmin, ymin, zmin);
  obb_points[1] = cp(xmax, ymin, zmin);
  obb_points[2] = cp(xmax, ymax, zmin);
  obb_points[3] = cp(xmin, ymax, zmin);

  obb_points[4] = cp(xmin, ymax, zmax); // see order in make_hexahedron()...
  obb_points[5] = cp(xmin, ymin, zmax);
  obb_points[6] = cp(xmax, ymin, zmax);
  obb_points[7] = cp(xmax, ymax, zmax);

  // Apply the inverse rotation to the rotated axis aligned bounding box
  for(std::size_t i=0; i<8; ++i)
  {
    obb_points[i] = inverse_transformation.transform(obb_points[i]);
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
    std::cout << "  OBB[" << i << "] = " << obb_points[i] << std::endl;
#endif
  }
}

template <typename PointRange, typename Traits>
void compute_best_transformation(const PointRange& points,
                                 typename Traits::Aff_transformation_3& transformation,
                                 typename Traits::Aff_transformation_3& inverse_transformation,
                                 CGAL::Random& rng,
                                 const Traits& traits)
{
  typedef typename Traits::Matrix                                    Matrix;
  typedef typename Traits::Aff_transformation_3                      Aff_transformation_3;

  CGAL_assertion(points.size() >= 3);

  const std::size_t max_generations = 100;
  const std::size_t population_size = 30;
  const std::size_t nelder_mead_iterations = 20;

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  CGAL::Real_timer timer;
  timer.start();
#endif

  Evolution<PointRange, Traits> search_solution(points, rng, traits);
  search_solution.evolve(max_generations, population_size, nelder_mead_iterations);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  std::cout << "evolve: " << timer.time() << std::endl;
  timer.reset();
#endif

  const Matrix& rot = search_solution.get_best_vertex().matrix();

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

// Following two functions are overloads to dispatch depending on return type
template <typename PointRange, typename K, typename Traits>
void construct_oriented_bounding_box(const PointRange& points,
                                     CGAL::Aff_transformation_3<K>& transformation,
                                     CGAL::Random& rng,
                                     const Traits& traits)
{
  typename Traits::Aff_transformation_3 inverse_transformation;
  compute_best_transformation(points, transformation, inverse_transformation, rng, traits);
}

template <typename PointRange, typename Array, typename Traits>
void construct_oriented_bounding_box(const PointRange& points,
                                     Array& obb_points,
                                     CGAL::Random& rng,
                                     const Traits& traits,
                                     typename boost::enable_if<
                                       typename boost::has_range_iterator<Array>
                                     >::type* = 0)
{
  typename Traits::Aff_transformation_3 transformation, inverse_transformation;
  compute_best_transformation(points, transformation, inverse_transformation, rng, traits);

  construct_oriented_bounding_box(points, transformation, inverse_transformation, obb_points, traits);
}

template <typename PointRange, typename PolygonMesh, typename Traits>
void construct_oriented_bounding_box(const PointRange& points,
                                     PolygonMesh& pm,
                                     CGAL::Random& rng,
                                     const Traits& traits,
                                     typename boost::disable_if<
                                       typename boost::has_range_iterator<PolygonMesh>
                                     >::type* = 0)
{
  typename Traits::Aff_transformation_3 transformation, inverse_transformation;
  compute_best_transformation(points, transformation, inverse_transformation, rng, traits);

  std::array<typename Traits::Point_3, 8> obb_points;
  construct_oriented_bounding_box(points, transformation, inverse_transformation, obb_points, traits);

  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], pm);
}

// Entry point, decide whether to compute the CH_3 or not
template <typename PointRange, typename Output, typename Traits>
void construct_oriented_bounding_box(const PointRange& points,
                                     const bool use_ch,
                                     Output& output,
                                     CGAL::Random& rng,
                                     const Traits& traits)
{
  typedef typename Traits::Point_3                                   Point;

  CGAL_static_assertion((std::is_same<typename boost::range_value<PointRange>::type, Point>::value));

  if(use_ch) // construct the convex hull to reduce the number of points
  {
    std::vector<Point> ch_points;

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  CGAL::Real_timer timer;
  timer.start();
#endif

    CGAL::Convex_hull_traits_3<Traits> CH_traits;
    extreme_points_3(points, std::back_inserter(ch_points), CH_traits);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
    std::cout << "CH time: " << timer.time() << std::endl;
#endif

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
    std::cout << ch_points.size() << " points on the convex hull" << std::endl;
#endif

    return construct_oriented_bounding_box(ch_points, output, rng, traits);
  }
  else
  {
    return construct_oriented_bounding_box(points, output, rng, traits);
  }
}

} // namespace internal
} // namespace Optimal_bounding_box

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
/// The input can be either a range of 3D points, or a polygon mesh.
///
/// \cgalHeading{Output}
///
/// The result of the algorithm can be retrieved as either:
/// - the best affine transformation \f${\mathcal R}_b\f$ that the algorithm has found;
/// - an array of eight points, representing the best oriented bounding box (\f${\mathcal B}_b\f$)
///   that the algorithm has constructed, which is related to \f$ {\mathcal R}_b\f$ as it is
///   the inverse transformation of the axis-aligned bounding box of the transformed point set.
///   The order of the points in the array is the same as in the function
///   \link PkgBGLHelperFct `CGAL::make_hexahedron()` \endlink,
///   which can be used to construct a mesh from these points.
/// - a model of `MutableFaceGraph`
///
/// Note that when returning an array of points, these points are constructed from the axis-aligned
/// bounding box and some precision loss should therefore be expected if a kernel not providing
/// exact constructions is used.
///
/// The algorithm is based on a paper by Chang, Gorissen, and Melchior \cgalCite{cgal:cgm-fobbo-11}.

/// \ingroup PkgOptimalBoundingBox_Oriented_bounding_box
///
/// The function `oriented_bounding_box` computes an approximation of the <i>optimal bounding box</i>,
/// which is defined as the rectangular box with smallest volume of all the rectangular boxes containing
/// the input points.
///
/// See \ref PkgOptimalBoundingBox_Oriented_bounding_box for more information.
///
/// \tparam PointRange a model of `Range`. The value type may not be equal to the type `%Point_3` of the traits class
///                    if a point map is provided via named parameters (see below) to access points.
/// \tparam Output either the type `Aff_transformation_3` of the traits class,
///                or `std::array<Point, 8>` with `Point` being equivalent to the type `%Point_3` of the traits class,
///                or a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param points the input range
/// \param out the resulting array of points or affine transformation
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{point_map}
///     \cgalParamDescription{a property map associating points to the elements of the point range}
///     \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
///     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a model of `OrientedBoundingBoxTraits_3`}
///     \cgalParamDefault{a default-constructed object of type `CGAL::Oriented_bounding_box_traits_3<K>`,
///                       where `K` is a kernel type deduced from the point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{use_convex_hull}
///     \cgalParamDescription{Parameter used in the construction of oriented bounding box to indicate
///                           whether the algorithm should first extract the extreme points (points
///                           that are on the 3D convex hull) of the input data range
///                           to accelerate the computation of the bounding box.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`true`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <typename PointRange,
          typename Output,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void oriented_bounding_box(const PointRange& points,
                           Output& out,
                           const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
                           , typename boost::enable_if<
                               typename boost::has_range_iterator<PointRange>
                           >::type* = 0
#endif
                           )
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type               PointMap;

#if defined(CGAL_EIGEN3_ENABLED)
  typedef typename boost::property_traits<PointMap>::value_type                         Point;
  typedef typename CGAL::Kernel_traits<Point>::type                                     K;
  typedef Oriented_bounding_box_traits_3<K>                                             Default_traits;
#else
  typedef void                                                                          Default_traits;
#endif

  typedef typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                       CGAL_BGL_NP_CLASS,
                                                       Default_traits>::type            Geom_traits;

  CGAL_static_assertion_msg(!(std::is_same<Geom_traits, void>::value),
                            "You must provide a traits class or have Eigen enabled!");

  Geom_traits traits = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  const bool use_ch = choose_parameter(get_parameter(np, internal_np::use_convex_hull), true);
  const unsigned int seed = choose_parameter(get_parameter(np, internal_np::random_seed), -1); // undocumented

  CGAL::Random fixed_seed_rng(seed);
  CGAL::Random& rng = (seed == unsigned(-1)) ? CGAL::get_default_random() : fixed_seed_rng;

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
  std::cout << "Random seed: " << rng.get_seed() << std::endl;
#endif

  // @todo handle those cases (or call min_rectangle_2 with a projection)
  if(points.size() <= 3)
  {
    std::cerr << "The oriented bounding box cannot (yet) be computed for a mesh with fewer than 4 vertices!\n";
    return;
  }

  return Optimal_bounding_box::internal::construct_oriented_bounding_box(
           CGAL::make_range(
             boost::make_transform_iterator(points.begin(), CGAL::Property_map_to_unary_function<PointMap>(point_map)),
             boost::make_transform_iterator(points.end(), CGAL::Property_map_to_unary_function<PointMap>(point_map))),
           use_ch, out, rng, traits);
}

/// \ingroup PkgOptimalBoundingBox_Oriented_bounding_box
///
/// Extracts the vertices of the mesh as a point range and calls the overload using points as input.
///
/// \tparam PolygonMesh a model of `VertexListGraph`
/// \tparam Output either the type `Aff_transformation_3` of the traits class,
///                or `std::array<Point, 8>` with `Point` being equivalent to the type `%Point_3` of the traits class,
///                or a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param pmesh the input mesh
/// \param out the resulting array of points or affine transformation
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     should be available for the vertices of `pmesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a model of `OrientedBoundingBoxTraits_3`}
///     \cgalParamDefault{a default-constructed object of type `CGAL::Oriented_bounding_box_traits_3<K>`,
///                       where `K` is a kernel type deduced from the point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{use_convex_hull}
///     \cgalParamDescription{Parameter used in the construction of oriented bounding box to indicate
///                           whether the algorithm should first extract the extreme points (points
///                           that are on the 3D convex hull) of the input data range
///                           to accelerate the computation of the bounding box.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`true`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <typename PolygonMesh,
          typename Output,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
void oriented_bounding_box(const PolygonMesh& pmesh,
                           Output& out,
                           const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
                           , typename boost::disable_if<
                              typename boost::has_range_iterator<PolygonMesh>
                           >::type* = 0
#endif
                           )
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, CGAL_BGL_NP_CLASS>::const_type  VPM;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  oriented_bounding_box(vertices(pmesh), out, np.point_map(vpm));
}

/// \cond SKIP_IN_MANUAL

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Convenience overloads
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Input /*range or mesh*/, typename Output /*transformation, array, or mesh*/>
void oriented_bounding_box(const Input& in, Output& out)
{
  return oriented_bounding_box(in, out, CGAL::parameters::all_default());
}

/// \endcond

} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_ORIENTED_BOUNDING_BOX_H
