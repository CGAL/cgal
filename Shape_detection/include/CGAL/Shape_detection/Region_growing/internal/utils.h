// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <set>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <typeinfo>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Named parameters.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Face graph includes.
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<typename GeomTraits>
  class Default_sqrt {

  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const {

      // TODO: This happens for circles and cylinders only! Maybe after my
      // precision cleaning in the new revision PR, this will be gone for all platforms.
      if (value < FT(0)) return FT(0); // clamp to zero
      const bool is_value_ok = (value >= FT(0));
      if (!is_value_ok) { // TODO: remove that!
        std::cout.precision(20);
        std::cout << "- wrong value: " << value << std::endl;
      }
      CGAL_precondition(is_value_ok);
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits,
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {

  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) {
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {

  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) {
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  struct Compare_scores {

    const std::vector<FT>& m_scores;
    Compare_scores(const std::vector<FT>& scores) :
    m_scores(scores)
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {
      CGAL_precondition(i < m_scores.size());
      CGAL_precondition(j < m_scores.size());
      return m_scores[i] > m_scores[j];
    }
  };

  template<
  typename Traits,
  typename InputRange,
  typename ItemMap>
  std::pair<typename Traits::Line_2, typename Traits::FT>
  create_line_2(
    const InputRange& input_range, const ItemMap item_map,
    const std::vector<std::size_t>& region, const Traits&) {

    using FT = typename Traits::FT;
    using Line_2 = typename Traits::Line_2;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_2 = typename ITraits::Point_2;
    using ILine_2 = typename ITraits::Line_2;

    using Input_type = typename ItemMap::value_type;
    using Item = typename std::conditional<
      std::is_same<typename Traits::Point_2, Input_type>::value,
      typename ITraits::Point_2, typename ITraits::Segment_2 >::type;

    std::vector<Item> items;
    CGAL_precondition(region.size() > 0);
    items.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const std::size_t item_index : region) {
      CGAL_precondition(item_index < input_range.size());
      const auto& key = *(input_range.begin() + item_index);
      const auto& item = get(item_map, key);
      items.push_back(iconverter(item));
    }
    CGAL_precondition(items.size() == region.size());

    ILine_2 fitted_line;
    IPoint_2 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_2(
      items.begin(), items.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<0>(), ITraits(),
      CGAL::Eigen_diagonalize_traits<IFT, 2>());

    const Line_2 line(
      static_cast<FT>(fitted_line.a()),
      static_cast<FT>(fitted_line.b()),
      static_cast<FT>(fitted_line.c()));
    return std::make_pair(line, static_cast<FT>(score));
  }

  template<
  typename Traits,
  typename InputRange,
  typename ItemMap>
  std::pair<typename Traits::Line_3, typename Traits::FT>
  create_line_3(
    const InputRange& input_range, const ItemMap item_map,
    const std::vector<std::size_t>& region, const Traits&) {

    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;
    using Point_3 = typename Traits::Point_3;
    using Direction_3 = typename Traits::Direction_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using ILine_3 = typename ITraits::Line_3;

    using Input_type = typename ItemMap::value_type;
    using Item = typename std::conditional<
      std::is_same<typename Traits::Point_3, Input_type>::value,
      typename ITraits::Point_3, typename ITraits::Segment_3 >::type;

    std::vector<Item> items;
    CGAL_precondition(region.size() > 0);
    items.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const std::size_t item_index : region) {
      CGAL_precondition(item_index < input_range.size());
      const auto& key = *(input_range.begin() + item_index);
      const auto& item = get(item_map, key);
      items.push_back(iconverter(item));
    }
    CGAL_precondition(items.size() == region.size());

    ILine_3 fitted_line;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      items.begin(), items.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<0>(), ITraits(),
      CGAL::Eigen_diagonalize_traits<IFT, 3>());

    const auto p = fitted_line.point(0);
    const auto d = fitted_line.direction();
    const Point_3 init(
      static_cast<FT>(p.x()),
      static_cast<FT>(p.y()),
      static_cast<FT>(p.z()));
    const Direction_3 direction(
      static_cast<FT>(d.dx()),
      static_cast<FT>(d.dy()),
      static_cast<FT>(d.dz()));
    const Line_3 line(init, direction);
    return std::make_pair(line, static_cast<FT>(score));
  }

  template<
  typename Traits,
  typename InputRange,
  typename ItemMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
  create_plane(
    const InputRange& input_range, const ItemMap item_map,
    const std::vector<std::size_t>& region, const Traits&) {

    using FT = typename Traits::FT;
    using Plane_3 = typename Traits::Plane_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using IPlane_3 = typename ITraits::Plane_3;

    using Input_type = typename ItemMap::value_type;
    using Item = typename std::conditional<
      std::is_same<typename Traits::Point_3, Input_type>::value,
      typename ITraits::Point_3, typename ITraits::Segment_3 >::type;

    std::vector<Item> items;
    CGAL_precondition(region.size() > 0);
    items.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const std::size_t item_index : region) {
      CGAL_precondition(item_index < input_range.size());
      const auto& key = *(input_range.begin() + item_index);
      const auto& item = get(item_map, key);
      items.push_back(iconverter(item));
    }
    CGAL_precondition(items.size() == region.size());

    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      items.begin(), items.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<0>(), ITraits(),
      CGAL::Eigen_diagonalize_traits<IFT, 3>());

    const Plane_3 plane(
      static_cast<FT>(fitted_plane.a()),
      static_cast<FT>(fitted_plane.b()),
      static_cast<FT>(fitted_plane.c()),
      static_cast<FT>(fitted_plane.d()));
    return std::make_pair(plane, static_cast<FT>(score));
  }

  template<
  typename Traits,
  typename FaceGraph,
  typename FaceRange,
  typename VertexToPointMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
  create_plane_from_faces(
    const FaceGraph& face_graph,
    const FaceRange& face_range,
    const VertexToPointMap vertex_to_point_map,
    const std::vector<std::size_t>& region, const Traits&) {

    using FT = typename Traits::FT;
    using Plane_3 = typename Traits::Plane_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using IPlane_3 = typename ITraits::Plane_3;

    std::vector<IPoint_3> points;
    CGAL_precondition(region.size() > 0);
    points.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const std::size_t face_index : region) {
      CGAL_precondition(face_index < face_range.size());
      const auto face = *(face_range.begin() + face_index);

      const auto hedge = halfedge(face, face_graph);
      const auto vertices = vertices_around_face(hedge, face_graph);
      CGAL_precondition(vertices.size() > 0);

      for (const auto vertex : vertices) {
        const auto& point = get(vertex_to_point_map, vertex);
        points.push_back(iconverter(point));
      }
    }
    CGAL_precondition(points.size() >= region.size());

    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      points.begin(), points.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<0>(), ITraits(),
      CGAL::Eigen_diagonalize_traits<IFT, 3>());

    const Plane_3 plane(
      static_cast<FT>(fitted_plane.a()),
      static_cast<FT>(fitted_plane.b()),
      static_cast<FT>(fitted_plane.c()),
      static_cast<FT>(fitted_plane.d()));
    return std::make_pair(plane, static_cast<FT>(score));
  }

  template<
  typename PointRange,
  typename Sqrt,
  typename Squared_distance_2,
  typename Point_2,
  typename FT>
  bool create_circle_2(
    const PointRange& points,
    const Sqrt& sqrt,
    const Squared_distance_2& squared_distance_2,
    Point_2& center, FT& radius) {

    using Diagonalize_traits = Eigen_diagonalize_traits<FT, 4>;
    using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;

    using Vector = typename Diagonalize_traits::Vector;
    using Matrix = typename Diagonalize_traits::Matrix;

    // Use bbox to compute diagonalization with smaller coordinates to
    // avoid loss of precision when inverting large coordinates.
    const Bbox_2 bbox = bbox_2(points.begin(), points.end());

    // Circle least squares fitting:
    // Circle of center (a, b) and radius R:
    // Ri = sqrt((xi - a)^2 + (yi - b)^2)
    // Minimize sum(Ri^2 - R^2)^2
    // -> Minimize sum(xi^2 + yi^2 − 2 * a * xi − 2 * b * yi + a^2 + b^2 − R^2)^2
    // Let B = -2a; C = -2b; D = a^2 + b^2 - R^2; and ri = x^2 + y^2.
    // -> Minimize sum(D + B * xi + C * yi + ri)^2
    // -> Minimize sum(1 + B / D * xi + C / D * yi + ri / D)^2
    // -> system of linear equations
    // -> diagonalize matrix
    //    NB   x   y   r
    //        xx  xy  xr
    //            yy  yr
    //                rr
    // -> center coordinates = [
    // -0.5 * eigenvector(1) / eigenvector(3);
    // -0.5 * eigenvector(2) / eigenvector(3); ]
    Covariance_matrix A = {
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0)
    };

    A[0] = static_cast<FT>(points.size());
    for (const Point_2& p : points) {

      const FT x = p.x() - bbox.xmin();
      const FT y = p.y() - bbox.ymin();
      const FT r = x * x + y * y;

      A[1] += x;
      A[2] += y;
      A[3] += r;
      A[4] += x * x;
      A[5] += x * y;
      A[6] += x * r;
      A[7] += y * y;
      A[8] += y * r;
      A[9] += r * r;
    }

    Vector eigenvalues = {
      FT(0), FT(0), FT(0), FT(0)
    };
    Matrix eigenvectors = {
      FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0)
    };
    Diagonalize_traits::
      diagonalize_selfadjoint_covariance_matrix(A, eigenvalues, eigenvectors);

    // Perfect line case, no circle can be fitted.
    if (eigenvectors[3] == 0) {
      return false;
    }

    // Other cases.
    const FT half = FT(1) / FT(2);
    center = Point_2(
      bbox.xmin() - half * (eigenvectors[1] / eigenvectors[3]),
      bbox.ymin() - half * (eigenvectors[2] / eigenvectors[3]));

    radius = FT(0);
    for (const Point_2& p : points) {
      radius += sqrt(squared_distance_2(p, center));
    }
    radius /= static_cast<FT>(points.size());
    return true;
  }

  template<
  typename PointRange,
  typename Sqrt,
  typename Squared_distance_3,
  typename Point_3,
  typename FT>
  bool create_sphere_3(
    const PointRange& points,
    const Sqrt& sqrt,
    const Squared_distance_3& squared_distance_3,
    Point_3& center, FT& radius) {

    using Diagonalize_traits = Eigen_diagonalize_traits<FT, 5>;
    using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;

    using Vector = typename Diagonalize_traits::Vector;
    using Matrix = typename Diagonalize_traits::Matrix;

    // Use bbox to compute diagonalization with smaller coordinates to
    // avoid loss of precision when inverting large coordinates.
    const Bbox_3 bbox = bbox_3(points.begin(), points.end());

    // Sphere least squares fitting.
    // See create_circle_2() above for more details on computation.
    Covariance_matrix A = {
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0)
    };

    A[0] = static_cast<FT>(points.size());
    for (const Point_3& p : points) {

      const FT x = p.x() - bbox.xmin();
      const FT y = p.y() - bbox.ymin();
      const FT z = p.z() - bbox.zmin();
      const FT r = x * x + y * y + z * z;

      A[1] += x;
      A[2] += y;
      A[3] += z;
      A[4] += r;
      A[5] += x * x;
      A[6] += x * y;
      A[7] += x * z;
      A[8] += x * r;
      A[9] += y * y;
      A[10] += y * z;
      A[11] += y * r;
      A[12] += z * z;
      A[13] += z * r;
      A[14] += r * r;
    }

    Vector eigenvalues = {
      FT(0), FT(0), FT(0), FT(0), FT(0)
    };
    Matrix eigenvectors = {
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0),
      FT(0), FT(0), FT(0), FT(0), FT(0)
    };
    Diagonalize_traits::
      diagonalize_selfadjoint_covariance_matrix(A, eigenvalues, eigenvectors);

    // Perfect plane case, no sphere can be fitted.
    if (eigenvectors[4] == 0) {
      return false;
    }

    // Other cases.
    const FT half = FT(1) / FT(2);
    center = Point_3(
      bbox.xmin() - half * (eigenvectors[1] / eigenvectors[4]),
      bbox.ymin() - half * (eigenvectors[2] / eigenvectors[4]),
      bbox.zmin() - half * (eigenvectors[3] / eigenvectors[4]));

    radius = FT(0);
    for (const Point_3& p : points) {
      radius += sqrt(squared_distance_3(p, center));
    }
    radius /= static_cast<FT>(points.size());
    return true;
  }

  template<
  typename PointRange,
  typename PointMap,
  typename NormalMap,
  typename Sqrt,
  typename Squared_distance_3,
  typename Line_3,
  typename FT>
  bool create_cylinder_3(
    const PointRange& points,
    PointMap point_map,
    NormalMap normal_map,
    const Sqrt& sqrt,
    const Squared_distance_3& squared_distance_3,
    Line_3& axis, FT& radius) {

    using Point_3 = typename boost::property_traits<PointMap>::value_type;
    using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

    const Point_3& ref = get(point_map, *(points.begin()));

    radius = FT(0);
    std::size_t nb = 0;
    Vector_3 mean_axis = CGAL::NULL_VECTOR;
    Point_3 point_on_axis = ORIGIN;

    for (std::size_t i = 0; i < points.size() - 1; ++i) {

      Vector_3 v0 = get(normal_map, *std::next(points.begin(), i));
      v0 = v0 / sqrt(v0 * v0);
      Vector_3 v1 = get(normal_map, *std::next(points.begin(), i + 1));
      v1 = v1 / sqrt(v1 * v1);
      Vector_3 axis = cross_product(v0, v1);
      if (sqrt(axis.squared_length()) < FT(1) / FT(100)) {
        continue;
      }
      axis = axis / sqrt(axis * axis);

      const Point_3& p0 = get(point_map, *std::next(points.begin(), i));
      const Point_3& p1 = get(point_map, *std::next(points.begin(), i + 1));

      Vector_3 xdir = v0 - axis * (v0 * axis);
      xdir = xdir / sqrt(xdir * xdir);

      Vector_3 ydir = cross_product(axis, xdir);
      ydir = ydir / sqrt (ydir * ydir);

      const FT v1x =  v1 * ydir;
      const FT v1y = -v1 * xdir;

      Vector_3 d(p0, p1);
      const FT ox = xdir * d;
      const FT oy = ydir * d;
      const FT ldist = v1x * ox + v1y * oy;

      FT r = ldist / v1x;
      Point_3 point = p0 + xdir * r;
      const Line_3 line(point, axis);
      point = line.projection(ref);
      point_on_axis = barycenter(point_on_axis, static_cast<FT>(nb), point, FT(1));
      r += abs(r);

      if (nb != 0 && axis * mean_axis < 0) {
        axis = -axis;
      }
      mean_axis = mean_axis + axis;
      ++nb;
    }

    if (nb == 0) {
      return false;
    }

    mean_axis = mean_axis / sqrt(mean_axis * mean_axis);
    axis = Line_3(point_on_axis, mean_axis);
    radius /= static_cast<FT>(nb);

    radius = FT(0);
    for (const auto& p : points) {
      const Point_3& p0 = get(point_map, p);
      radius += sqrt(squared_distance_3(p0, axis));
    }
    radius /= static_cast<FT>(points.size());
    return true;
  }

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
