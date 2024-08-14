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
#include <CGAL/Shape_detection/Region_growing/internal/cylinder_fitting.h>

// STL includes.
#include <set>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <typeinfo>

// Boost headers.
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Named parameters.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iterator_range.h>
#ifdef CGAL_SD_RG_USE_PMP
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#else
#include <CGAL/centroid.h>
#endif

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<class T>
  struct reference_iterator_generator {
    using result_type = T;
    reference_iterator_generator(T iterator) : it(iterator) {  }
    result_type operator() () const {
      return it++;
    }
    mutable T it;
  };

  // TODO: this should be customizable in named function parameters
  template<class T, bool = CGAL::is_iterator<T>::value>
  struct hash_item {};

  template<class T>
  struct hash_item<T, false> {
    std::size_t operator()(T i) const {
      using boost::hash_value;
      return hash_value(i);
    }
  };

  template<class T>
  struct hash_item<T, true> {
    std::size_t operator()(T i) const {
      using boost::hash_value;
      return boost::hash_value(i.operator->());
    }
  };

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
  typename Region,
  typename PrimitiveMap>
  std::pair<typename Traits::Line_2, typename Traits::FT>
  create_line_2(
    const Region& region, const PrimitiveMap primitive_map, const Traits&) {

    using FT = typename Traits::FT;
    using Line_2 = typename Traits::Line_2;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = CGAL::Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_2 = typename ITraits::Point_2;
    using ILine_2 = typename ITraits::Line_2;

    using Item = typename boost::property_traits<PrimitiveMap>::key_type;
    using Primitive = typename boost::property_traits<PrimitiveMap>::value_type;
    using EPIC_Primitive = typename std::conditional<
      std::is_same<typename Traits::Point_2, Primitive>::value,
      typename ITraits::Point_2, typename ITraits::Segment_2 >::type;

    std::vector<EPIC_Primitive> elements;
    CGAL_precondition(region.size() > 0);
    elements.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (Item item : region) {
      const Primitive& element = get(primitive_map, item);
      elements.push_back(iconverter(element));
    }
    CGAL_precondition(elements.size() == region.size());

    ILine_2 fitted_line;
    IPoint_2 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_2(
      elements.begin(), elements.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<Feature_dimension<Primitive>::value>(), ITraits(),
      CGAL::Eigen_diagonalize_traits<IFT, 2>());

    const Line_2 line(
      static_cast<FT>(fitted_line.a()),
      static_cast<FT>(fitted_line.b()),
      static_cast<FT>(fitted_line.c()));

    return std::make_pair(line, static_cast<FT>(score));
  }

  template<
  typename Traits,
  typename Region,
  typename PrimitiveMap>
  std::pair<typename Traits::Line_3, typename Traits::FT>
  create_line_3(
    const Region& region, const PrimitiveMap primitive_map, const Traits&) {

    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;
    using Point_3 = typename Traits::Point_3;
    using Direction_3 = typename Traits::Direction_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = CGAL::Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using ILine_3 = typename ITraits::Line_3;

    using Item = typename boost::property_traits<PrimitiveMap>::key_type;
    using Primitive = typename boost::property_traits<PrimitiveMap>::value_type;
    using EPIC_Primitive = typename std::conditional<
      std::is_same<typename Traits::Point_3, Primitive>::value,
      typename ITraits::Point_3, typename ITraits::Segment_3 >::type;

    std::vector<EPIC_Primitive> elements;
    CGAL_precondition(region.size() > 0);
    elements.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (Item item : region) {
      const Primitive& element = get(primitive_map, item);
      elements.push_back(iconverter(element));
    }
    CGAL_precondition(elements.size() == region.size());

    ILine_3 fitted_line;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      elements.begin(), elements.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<Feature_dimension<Primitive>::value>(), ITraits(),
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
  typename Region,
  typename ItemMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
  create_plane(
    const Region& region, const ItemMap item_map, const Traits&) {

    using FT = typename Traits::FT;
    using Plane_3 = typename Traits::Plane_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = CGAL::Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using IPlane_3 = typename ITraits::Plane_3;

    using Input_type = typename ItemMap::value_type;
    using Element = typename std::conditional<
      std::is_same<typename Traits::Point_3, Input_type>::value,
      typename ITraits::Point_3, typename ITraits::Segment_3 >::type;

    std::vector<Element> elements;
    CGAL_precondition(region.size() > 0);
    elements.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (auto item : region) {
      const auto& element = get(item_map, item);
      elements.push_back(iconverter(element));
    }
    CGAL_precondition(elements.size() == region.size());

    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      elements.begin(), elements.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<Feature_dimension<Element>::value>(), ITraits(),
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
    typename InputRange,
    typename Region,
    typename ItemMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
    create_plane(
      const InputRange &input_range,
      const ItemMap item_map, const Region& region, const Traits& traits) {
    std::vector<typename InputRange::const_iterator> tmp;
    tmp.resize(region.size());

    for (std::size_t idx : region)
      tmp[idx] = input_range.begin() + idx;

    return create_plane(tmp, item_map, traits);
  }

  template <class Traits>
  void
  triangulate_face(const std::vector<typename Traits::Point_3>& points,
                         std::vector<typename Traits::Triangle_3>& triangles)
  {
#ifdef CGAL_SD_RG_USE_PMP
    std::vector<CGAL::Triple<int, int, int>> output;

    Polygon_mesh_processing::triangulate_hole_polyline(points, std::back_inserter(output), parameters::use_2d_constrained_delaunay_triangulation(true));

    triangles.reserve(output.size());
    for (const auto& t : output)
      triangles.emplace_back(points[t.first], points[t.second], points[t.third]);
#else
    //use a triangulation using the centroid
    std::size_t nb_edges = points.size();
    typename Traits::Point_3 c = CGAL::centroid(points.begin(), points.end());
    triangles.reserve(nb_edges);
    for (std::size_t i=0; i<nb_edges-1; ++i)
      triangles.emplace_back(points[i], points[i+1], c);
    triangles.emplace_back(points.back(), points.front(), c);
#endif
  }

  template<
  typename Traits,
  typename FaceGraph,
  typename Region,
  typename VertexToPointMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
  create_plane_from_faces(
    const FaceGraph& face_graph,
    const Region& region,
    const VertexToPointMap vertex_to_point_map, const Traits&) {

    using FT = typename Traits::FT;
    using Plane_3 = typename Traits::Plane_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = CGAL::Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using ITriangle_3 = typename ITraits::Triangle_3;
    using IPlane_3 = typename ITraits::Plane_3;

    std::vector<ITriangle_3> triangles;
    CGAL_precondition(region.size() > 0);
    triangles.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const typename Region::value_type face : region) {
      const auto hedge = halfedge(face, face_graph);
      const auto vertices = vertices_around_face(hedge, face_graph);
      CGAL_precondition(vertices.size() > 0);

      std::vector<IPoint_3> points;
      for (const auto vertex : vertices) {
        const auto& point = get(vertex_to_point_map, vertex);
        points.push_back(iconverter(point));
      }
      if (points.size()==3)
        triangles.push_back(ITriangle_3(points[0], points[1], points[2]));
      else
        triangulate_face<ITraits>(points, triangles);
    }
    CGAL_precondition(triangles.size() >= region.size());
    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      triangles.begin(), triangles.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<2>(), ITraits(),
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
    typename Region,
    typename FaceToTrianglesMap>
  std::pair<typename Traits::Plane_3, typename Traits::FT>
    create_plane_from_triangulated_faces(
      const Region& region,
      const FaceToTrianglesMap &face_to_triangles_map, const Traits&) {

    using FT = typename Traits::FT;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_3 = typename Traits::Triangle_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IConverter = CGAL::Cartesian_converter<Traits, ITraits>;

    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using ITriangle_3 = typename ITraits::Triangle_3;
    using IPlane_3 = typename ITraits::Plane_3;

    std::vector<ITriangle_3> triangles;
    CGAL_precondition(region.size() > 0);
    triangles.reserve(region.size());
    const IConverter iconverter = IConverter();

    for (const typename Region::value_type face : region) {
      const std::vector<Triangle_3>& tris = get(face_to_triangles_map, face);

      // Degenerate polygons are omitted.
      if (tris.empty())
        continue;

      for (const auto &tri : tris)
        triangles.push_back(iconverter(tri));
    }
    CGAL_precondition(!triangles.empty());
    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    const IFT score = CGAL::linear_least_squares_fitting_3(
      triangles.begin(), triangles.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<2>(), ITraits(),
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
  typename Region,
  typename PointMap>
  std::pair<std::pair<typename Traits::FT, typename Traits::Point_2>,
    typename Traits::FT>
  create_circle_2(
    const Region &region, const PointMap point_map,
    const Traits&, const bool compute_score) {

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    typename Get_sqrt<Traits>::Sqrt sqrt;
    typename Traits::Compute_squared_distance_2 squared_distance_2;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IFT = typename ITraits::FT;

    using Diagonalize_traits = Eigen_diagonalize_traits<IFT, 4>;
    using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;

    using Vector = typename Diagonalize_traits::Vector;
    using Matrix = typename Diagonalize_traits::Matrix;

    // Use bbox to compute diagonalization with smaller coordinates to
    // avoid loss of precision when inverting large coordinates.
    FT xmin = +FT(1000000000000);
    FT ymin = +FT(1000000000000);
    for (auto item : region) {
      const auto& point = get(point_map, item);
      xmin = (CGAL::min)(xmin, point.x());
      ymin = (CGAL::min)(ymin, point.y());
    }

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
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0)
    };

    A[0] = static_cast<IFT>(region.size());
    for (auto &item : region) {
      const auto& point = get(point_map, item);

      const IFT x = static_cast<IFT>(CGAL::to_double(point.x() - xmin));
      const IFT y = static_cast<IFT>(CGAL::to_double(point.y() - ymin));
      const IFT r = x * x + y * y;

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
      IFT(0), IFT(0), IFT(0), IFT(0)
    };
    Matrix eigenvectors = {
      IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0)
    };
    Diagonalize_traits::
      diagonalize_selfadjoint_covariance_matrix(A, eigenvalues, eigenvectors);

    // Perfect line case, no circle can be fitted.
    if (eigenvectors[3] == IFT(0)) {
      return std::make_pair(
        std::make_pair(FT(-1), Point_2()),
        -static_cast<FT>(((std::numeric_limits<double>::max))()));
    }
    CGAL_assertion(eigenvectors[3] != IFT(0));

    // Other cases.
    const FT half = FT(1) / FT(2);
    const Point_2 fitted_center = Point_2(
      xmin - half * static_cast<FT>(eigenvectors[1] / eigenvectors[3]),
      ymin - half * static_cast<FT>(eigenvectors[2] / eigenvectors[3])
    );

    FT fitted_radius = FT(0);
    for (auto &item : region) {
      const auto& point = get(point_map, item);
      fitted_radius += sqrt(squared_distance_2(point, fitted_center));
    }
    fitted_radius /= static_cast<FT>(region.size());

    // Compute score.
    FT score = FT(-1);
    if (compute_score) {
      score = FT(0);
      for (auto &item : region) {
        const auto& point = get(point_map, item);
        score -= CGAL::abs(sqrt(squared_distance_2(point, fitted_center))
                           - fitted_radius);
      }
    }

    return std::make_pair(
      std::make_pair(fitted_radius, fitted_center), score);
  }

  template<
  typename Traits,
  typename Region,
  typename PointMap>
  std::pair<std::pair<typename Traits::FT, typename Traits::Point_3>,
            typename Traits::FT>
  create_sphere(
    const Region &region, const PointMap point_map,
    const Traits&, const bool compute_score) {

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    typename Get_sqrt<Traits>::Sqrt sqrt;
    typename Traits::Compute_squared_distance_3 squared_distance_3;

    using ITraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using IFT = typename ITraits::FT;

    using Diagonalize_traits = Eigen_diagonalize_traits<IFT, 5>;
    using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;

    using Vector = typename Diagonalize_traits::Vector;
    using Matrix = typename Diagonalize_traits::Matrix;

    // Use bbox to compute diagonalization with smaller coordinates to
    // avoid loss of precision when inverting large coordinates.
    FT xmin = +FT(1000000000000);
    FT ymin = +FT(1000000000000);
    FT zmin = +FT(1000000000000);
    for (auto item : region) {
      const auto& point = get(point_map, item);
      xmin = (CGAL::min)(xmin, point.x());
      ymin = (CGAL::min)(ymin, point.y());
      zmin = (CGAL::min)(zmin, point.z());
    }

    // Sphere least squares fitting.
    // See create_circle_2() above for more details on computation.
    Covariance_matrix A = {
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0)
    };

    A[0] = static_cast<IFT>(region.size());
    for (auto item : region) {
      const auto& point = get(point_map, item);

      const IFT x = static_cast<IFT>(CGAL::to_double(point.x() - xmin));
      const IFT y = static_cast<IFT>(CGAL::to_double(point.y() - ymin));
      const IFT z = static_cast<IFT>(CGAL::to_double(point.z() - zmin));
      const IFT r = x * x + y * y + z * z;

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
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0)
    };
    Matrix eigenvectors = {
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0),
      IFT(0), IFT(0), IFT(0), IFT(0), IFT(0)
    };
    Diagonalize_traits::
      diagonalize_selfadjoint_covariance_matrix(A, eigenvalues, eigenvectors);

    // Perfect plane case, no sphere can be fitted.
    if (eigenvectors[4] == IFT(0)) {
      return std::make_pair(
        std::make_pair(FT(-1), Point_3()),
        -static_cast<FT>((std::numeric_limits<double>::max)()));
    }
    CGAL_assertion(eigenvectors[4] != IFT(0));

    // Other cases.
    const FT half = FT(1) / FT(2);
    const Point_3 fitted_center = Point_3(
      xmin - half * static_cast<FT>(eigenvectors[1] / eigenvectors[4]),
      ymin - half * static_cast<FT>(eigenvectors[2] / eigenvectors[4]),
      zmin - half * static_cast<FT>(eigenvectors[3] / eigenvectors[4])
    );

    FT fitted_radius = FT(0);
    for (auto item : region) {
      const auto& point = get(point_map, item);
      fitted_radius += sqrt(squared_distance_3(point, fitted_center));
    }
    fitted_radius /= static_cast<FT>(region.size());

    // Compute score.
    FT score = FT(-1);
    if (compute_score) {
      score = FT(0);
      for (auto item : region) {
        const auto& point = get(point_map, item);
        score -= CGAL::abs(sqrt(squared_distance_3(point, fitted_center))
                           - fitted_radius);
      }
    }

    return std::make_pair(
      std::make_pair(fitted_radius, fitted_center), score);
  }

  template<
  typename Traits,
  typename Region,
  typename PointMap,
  typename NormalMap>
  std::pair<std::pair<typename Traits::FT, typename Traits::Line_3>,
            typename Traits::FT>
    create_cylinder(
      const Region& region, const PointMap point_map,
      const NormalMap normal_map,
      const Traits& traits) {

    if (region.size() < 6)
      return std::make_pair(
        std::make_pair<typename Traits::FT, typename Traits::Line_3>
        (-1.0, typename Traits::Line_3()),
        (std::numeric_limits<double>::max)());

    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;

    Line_3 fitted_axis;
    FT squared_radius;

    FT error = fit_cylinder<Traits, Region, PointMap, NormalMap>
      (region, point_map, normal_map, fitted_axis,
        squared_radius, traits);

    // A negative squared_radius is returned if the cylinder fitting failed.
    // This can be the case if the normals are very close to each other.
    if (squared_radius < 0)
      return std::make_pair(
        std::make_pair<typename Traits::FT, typename Traits::Line_3>
        (-1.0, typename Traits::Line_3()),
        (std::numeric_limits<double>::max)());

    return std::make_pair(
      std::make_pair(sqrt(squared_radius), fitted_axis), error);
  }

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
