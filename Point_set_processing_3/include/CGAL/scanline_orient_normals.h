// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SCANLINE_ORIENT_NORMALS_H
#define CGAL_SCANLINE_ORIENT_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <boost/iterator/transform_iterator.hpp>

#if defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

//#define CGAL_SCANLINE_ORIENT_VERBOSE

namespace CGAL
{

/// \cond SKIP_IN_MANUAL
namespace Point_set_processing_3
{

namespace internal
{

template <typename Vector_3>
const Vector_3& vertical_vector()
{
  static Vector_3 v(0, 0, 1);
  return v;
}

template <typename Iterator, typename PointMap,
          typename ScanlineIDMap>
bool is_end_of_scanline (Iterator scanline_begin, Iterator it,
                         PointMap,
                         ScanlineIDMap scanline_id_map,
                         const Tag_false&) // no fallback
{
  return (get (scanline_id_map, *scanline_begin)
          != get (scanline_id_map, *it));
}

template <typename Iterator, typename PointMap,
          typename ScanlineIDMap>
bool is_end_of_scanline (Iterator scanline_begin, Iterator it,
                         PointMap point_map,
                         ScanlineIDMap,
                         const Tag_true&) // fallback
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename Kernel_traits<Point_3>::Kernel::Vector_3;

  if (std::distance (scanline_begin, it) < 3)
    return false;

  Iterator n_minus_1 = it; -- n_minus_1;
  Iterator n_minus_2 = n_minus_1; -- n_minus_2;
  Iterator n_minus_3 = n_minus_2; -- n_minus_3;

  const Point_3& p_minus_1 = get(point_map, *n_minus_1);
  const Point_3& p_minus_2 = get(point_map, *n_minus_2);
  const Point_3& p_minus_3 = get(point_map, *n_minus_3);

  // End of scanline reached if inversion of direction
  Vector_3 v32 (p_minus_3, p_minus_2);
  v32 = Vector_3 (v32.x(), v32.y(), 0);
  Vector_3 v21 (p_minus_2, p_minus_1);
  v21 = Vector_3 (v21.x(), v21.y(), 0);

  return (v32 * v21 < 0);
}

template <typename Iterator, typename PointMap>
std::pair<typename Kernel_traits<typename boost::property_traits
                                 <PointMap>::value_type>::Kernel::Vector_3,
          typename Kernel_traits<typename boost::property_traits
                                 <PointMap>::value_type>::Kernel::Vector_3>
scanline_base (Iterator begin, Iterator end,
               PointMap point_map)
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Kernel = typename Kernel_traits<Point_3>::Kernel;
  using Vector_3 = typename Kernel::Vector_3;

  using Line_3 = typename Kernel::Line_3;
  using Plane_3 = typename Kernel::Plane_3;

  const double limit = CGAL_PI * 30. / 180.;

  Line_3 pca2;
  linear_least_squares_fitting_3
    (boost::make_transform_iterator
     (begin, Property_map_to_unary_function<PointMap>(point_map)),
      boost::make_transform_iterator
     (end, Property_map_to_unary_function<PointMap>(point_map)),
      pca2, Dimension_tag<0>());

  Vector_3 pca_direction = pca2.to_vector();
  pca_direction = Vector_3 (pca_direction.x(), pca_direction.y(), 0);
  pca_direction = pca_direction / CGAL::approximate_sqrt(pca_direction.squared_length());

  Plane_3 pca3;
  linear_least_squares_fitting_3
    (boost::make_transform_iterator
     (begin, Property_map_to_unary_function<PointMap>(point_map)),
     boost::make_transform_iterator
     (end, Property_map_to_unary_function<PointMap>(point_map)),
     pca3, Dimension_tag<0>());

  Vector_3 orthogonal = pca3.orthogonal_vector();
  Vector_3 vertical = CGAL::cross_product (pca_direction, orthogonal);
  if (vertical * vertical_vector<Vector_3>() < 0)
    vertical = -vertical;

  vertical = vertical / CGAL::approximate_sqrt(vertical.squared_length());

  if (std::acos(vertical * vertical_vector<Vector_3>()) < limit)
    return std::make_pair (pca_direction, vertical);

  // if plane diverges from the vertical more than 30 degrees, then
  // fallback to 0 0 1 vertical vector

  // Dummy begin->end vector version
  Iterator last = end; -- last;
  const Point_3& pbegin = get (point_map, *begin);
  const Point_3& plast = get (point_map, *last);
  Vector_3 direction (Point_3 (pbegin.x(), pbegin.y(), 0),
                      Point_3 (plast.x(), plast.y(), 0));
  direction = direction / CGAL::approximate_sqrt (direction.squared_length());

  return std::make_pair (direction, vertical_vector<Vector_3>());
}

template <typename Vector_3>
bool normal_along_scanline_is_inverted
(const Vector_3& normal, const Vector_3& line_of_sight)
{
  return (line_of_sight * normal < 0);
}

template <typename Iterator, typename PointMap, typename Vector_3>
std::pair<typename boost::property_traits<PointMap>::value_type, bool>
estimate_scan_position (Iterator begin, Iterator end, PointMap point_map,
                        const std::vector<Vector_3>& lines_of_sight)
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  typedef Eigen::Matrix3d Matrix;
  typedef Eigen::Vector3d Vector;
  typedef Eigen::ConjugateGradient<Matrix> Solver;

  Matrix R;
  R << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  Vector Q;
  Q << 0, 0, 0;

  std::size_t idx = 0;
  for (Iterator it = begin; it != end; ++ it, ++ idx)
  {
    const Point_3& p = get (point_map, *it);
    const Vector_3& v = lines_of_sight[idx];

    Vector n;
    n << v.x(), v.y(), v.z();

    Matrix I_nnt = Matrix::Identity() - n * n.transpose();
    Vector a;
    a << p.x(), p.y(), p.z();

    R += I_nnt;
    Q += I_nnt * a;
  }

  Solver solver(R);

  if (solver.info() != Eigen::Success)
    return std::make_pair (ORIGIN, false);

  Vector p = solver.solve(Q);
  if (solver.info() != Eigen::Success)
    return std::make_pair (ORIGIN, false);

  return std::make_pair (Point_3(p(0), p(1), p(2)), true);
}


template <typename Iterator, typename PointMap, typename NormalMap,
          typename ScanAngleMap>
void orient_scanline (Iterator begin, Iterator end,
                      PointMap point_map,
                      NormalMap normal_map,
                      ScanAngleMap scan_angle_map,
                      const Tag_false&) // no fallback scan angle
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  Vector_3 direction;
  Vector_3 vertical;
  std::tie (direction, vertical)
    = scanline_base (begin, end, point_map);

  std::vector<Vector_3> lines_of_sight;
  lines_of_sight.reserve (std::distance (begin, end));

  double mean_z = 0;
  for (Iterator it = begin; it != end; ++ it)
  {
    double angle = CGAL_PI * static_cast<double>(get (scan_angle_map, *it)) / 180.;
    Vector_3 los = direction * std::sin(angle) + vertical * std::cos(angle);
    lines_of_sight.push_back (los);
    mean_z += get (point_map, *it).z();
  }
  mean_z /= std::distance (begin, end);

#ifdef CGAL_SCANORIENT_DUMP_RANDOM_SCANLINES
  if (rand() % 1000 == 0 && std::distance(begin, end) > 10)
  {
    std::ofstream ofile ("scanline.polylines.txt");
    ofile.precision(18);
    ofile << std::distance (begin, end);
    for (Iterator it = begin; it != end; ++ it)
      ofile << " " << get(point_map, *it);
    ofile << std::endl;

    std::ofstream ofile2 ("base.polylines.txt");
    Point_3 orig = get (point_map, *(begin + std::distance(begin, end) / 2));
    double dist = CGAL::approximate_sqrt (CGAL::squared_distance
                                          (get(point_map, *begin), orig));

    ofile2.precision(18);
    ofile2 << "2 " << orig << " " << orig + direction * dist << std::endl
           << "2 " << orig << " " << orig + vertical * dist << std::endl;
  }
#endif

  Point_3 scan_position;
  bool solver_success;
  std::tie (scan_position, solver_success)
    = estimate_scan_position (begin, end, point_map, lines_of_sight);

  // If solver failed OR if scan position is detected under the
  // scanline (obviously wrong)
  if (!solver_success || scan_position.z() < mean_z)
  {
#ifdef CGAL_SCANLINE_ORIENT_VERBOSE
    if (!solver_success)
      std::cerr << "Inverting because olver failed: ";
    else
      std::cerr << "Inverting because scanner under scanline: ";
#endif

    direction = -direction;
    std::size_t idx = 0;
    for (Iterator it = begin; it != end; ++ it, ++ idx)
    {
      double angle = CGAL_PI * static_cast<double>(get (scan_angle_map, *it)) / 180.;
      Vector_3 los = direction * std::sin(angle) + vertical * std::cos(angle);
      lines_of_sight[idx] = los;
    }

    std::tie (scan_position, solver_success)
      = estimate_scan_position (begin, end, point_map, lines_of_sight);

#ifdef CGAL_SCANLINE_ORIENT_VERBOSE
    if (solver_success && scan_position.z() > mean_z)
      std::cerr << "SOLVED" << std::endl;
    else if (!solver_success)
      std::cerr << "FAILED, solver failure" << std::endl;
    else
      std::cerr << "FAILED, scanner under scanline" << std::endl;
#endif
  }

  std::size_t idx = 0;
  for (Iterator it = begin; it != end; ++ it, ++ idx)
  {
    const Vector_3 los = lines_of_sight[idx];
    const Vector_3& normal = get (normal_map, *it);
    if (normal_along_scanline_is_inverted (normal, los))
      put (normal_map, *it, -normal);
  }
}

template <typename Iterator, typename PointMap, typename NormalMap,
          typename ScanAngleMap>
void orient_scanline (Iterator begin, Iterator end,
                      PointMap point_map,
                      NormalMap normal_map,
                      ScanAngleMap,
                      const Tag_true&) // fallback scan angle
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  Vector_3 direction;
  Vector_3 vertical;
  std::tie (direction, vertical)
    = scanline_base (begin, end, point_map);

  // Estimate scanner position:
  // average XY-projected point, located above
  double mean_x = 0.;
  double mean_y = 0.;
  double max_z = -(std::numeric_limits<double>::max)();
  std::size_t nb = 0;

  for (Iterator it = begin; it != end; ++ it)
  {
    const Point_3& p = get (point_map, *it);
    mean_x += p.x();
    mean_y += p.y();
    max_z = (std::max)(max_z, p.z());
    ++ nb;
  }

  Iterator last = end; -- last;
  double length = CGAL::approximate_sqrt (CGAL::squared_distance
                                          (get (point_map, *begin),
                                           get (point_map, *last)));

  Point_3 scan_position (mean_x / nb, mean_y / nb,
                         max_z + length * 2);

  for (Iterator it = begin; it != end; ++ it)
  {
    Vector_3 line_of_sight (get(point_map, *it), scan_position);
    const Vector_3& normal = get (normal_map, *it);
    if (normal_along_scanline_is_inverted (normal, line_of_sight))
      put (normal_map, *it, -normal);
  }
}

} // namespace internal

} // namespace Point_set_processing_3
/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms

   orients the normals of the range of `points` by estimating a line
   of sight and checking its consistency with the current normal orientation.

   \warning This function requires the input `points` to be ordered
   along scanlines aligned on the XY-plane. It is typically designed
   for 2.5D urban datasets acquired through, for example, airborne
   LIDAR devices.

   First, scanlines are estimated as subranges of `points` by
   iterating on `points`:

   - if the named parameter `scanline_id_map` is provided, the range
     is cutted everytime the id changes.

   - if no scanline ID map is provided, a fallback method simply cuts
     the range everytime 3 consecutive points form an acute angle on
     the projected XY-plane. This fallback method gives suboptimal
     results.

   Then, the line of sight (estimated vector between a point and the
   position of the scanner at its time of acquisition) is estimated:

   - if `scan_angle` is provided, the line of sight can be directly
     computed as a combination of the estimated scanline and of the
     scan angle.

   - if no scan angle map is provided, then for each scanline, the
     position of the scanner is estimated as being above of the
     barycenter of the points of the scanline projected on the
     XY-plane. This fallback method gives suboptimal results.

   Once the line of sight is estimated for each point, the normals are
   oriented by checking, for each of them, if the line of sight and the
   normal vector give a positive scalar product. If they don't, then
   the normal vector is inverted.

   \note This method gives optimal results when `scanline_id_map`
   and `scan_angle` are provided. Correct results may still be
   produced in the absence of either one or both of these properties,
   as long as the point set is ordered in 2.5D scanlines.

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the
       elements of the point set `points`}
       \cgalParamType{a model of `WritablePropertyMap` whose key type
                      is the value type of the iterator of
                      `PointRange` and whose value type is
                      `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{scan_angle_map}
       \cgalParamDescription{a property map associating the angle of
       acquisition (in degrees) to the elements of the point set
       `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type
                      is the value type of the iterator of
                      `PointRange` and whose value type is convertible
                      to `double`}
     \cgalParamNEnd

     \cgalParamNBegin{scanline_id_map}
       \cgalParamDescription{a property map associating a scanline ID
       to the elements of the point set `points`. A scanline is
       detected as a consecutive subrange of items in the input range
       `point` whose ID are identical. IDs do not need to be unique,
       they just need to be different for two consecutive
       scanlines. The LAS property `scan_direction_flag` (whose values
       are either 0 or 1 depending on the direction of the scanner)
       can be used.}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type
                      is the value type of the iterator of
                      `PointRange` and whose value type is a model of
                      `EqualityComparable`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd
*/
template <typename PointRange, typename NamedParameters>
void scanline_orient_normals (PointRange& points, const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Iterator = typename PointRange::iterator;
  using Value_type = typename std::iterator_traits<Iterator>::value_type;
  using PointMap = typename CGAL::GetPointMap<PointRange, NamedParameters>::type;
  using NormalMap = typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type;

  using No_map = Constant_property_map<Value_type, int>;

  using ScanAngleMap = typename internal_np::Lookup_named_param_def
    <internal_np::scan_angle_t, NamedParameters, No_map>::type;
  using Fallback_scan_angle = Boolean_tag<std::is_same<ScanAngleMap, No_map>::value>;

  using ScanlineIDMap = typename internal_np::Lookup_named_param_def
    <internal_np::scanline_id_t, NamedParameters, No_map>::type;
  using Fallback_scanline_ID = Boolean_tag<std::is_same<ScanlineIDMap, No_map>::value>;

  CGAL_static_assertion_msg(!(std::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap
                              <PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));
  ScanAngleMap scan_angle_map = choose_parameter<ScanAngleMap>
    (get_parameter(np, internal_np::scan_angle_map));
  ScanlineIDMap scanline_id_map = choose_parameter<ScanlineIDMap>
    (get_parameter(np, internal_np::scanline_id_map));

  std::size_t nb_scanlines = 1;

#ifdef CGAL_SCANORIENT_DUMP_RANDOM_SCANLINES
  std::ofstream ofile ("scanlines.polylines.txt");
  ofile.precision(18);
#endif

  CGAL::Bbox_3 bbox = CGAL::bbox_3
    (boost::make_transform_iterator
     (points.begin(), Property_map_to_unary_function<PointMap>(point_map)),
     boost::make_transform_iterator
     (points.end(), Property_map_to_unary_function<PointMap>(point_map)));

  double bbox_diagonal
    = CGAL::approximate_sqrt((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
                             + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
                             + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));
  double limit = 0.05 * bbox_diagonal;

  Iterator scanline_begin = points.begin();
  for (Iterator it = points.begin(); it != points.end(); ++ it)
  {
    bool force_cut = false;
    if (it != points.begin())
    {
      Iterator prev = it; -- prev;
      force_cut = (CGAL::squared_distance (get (point_map, *prev),
                                           get (point_map, *it)) > limit * limit);
    }

    if (Point_set_processing_3::internal::is_end_of_scanline
        (scanline_begin, it, point_map, scanline_id_map,
         Fallback_scanline_ID()) || force_cut)
    {
      Point_set_processing_3::internal::orient_scanline
        (scanline_begin, it, point_map, normal_map,
         scan_angle_map, Fallback_scan_angle());

#ifdef CGAL_SCANORIENT_DUMP_RANDOM_SCANLINES
      ofile << std::distance (scanline_begin, it);
      for (Iterator it2 = scanline_begin; it2 != it; ++ it2)
        ofile << " " << get (point_map, *it2);
      ofile << std::endl;
#endif

      scanline_begin = it;
      ++ nb_scanlines;
    }
  }

  Point_set_processing_3::internal::orient_scanline
    (scanline_begin, points.end(), point_map, normal_map,
     scan_angle_map, Fallback_scan_angle());

#ifdef CGAL_SCANLINE_ORIENT_VERBOSE
  std::cerr << nb_scanlines << " scanline(s) identified (mean length = "
            << std::size_t(points.size() / double(nb_scanlines))
            << " point(s))" << std::endl;
#endif
}

template <typename PointRange>
void scanline_orient_normals (PointRange& points)
{
  return scanline_orient_normals (points,
                               CGAL::Point_set_processing_3::parameters::all_default(points));
}

} // namespace CGAL

#endif // CGAL_EIGEN3_ENABLED

#endif // CGAL_SCANLINE_ORIENT_NORMALS_H
