// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_OPTIMIZE_2_H
#define CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_OPTIMIZE_2_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_type_config.h>

#include <iostream>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

enum PROJECTION_DIRECTION
{
  ALONG_X = 0,
  ALONG_Y,
  ALONG_Z
};

// Now, we would like to do all of this with projection traits... Unfortunately, it's missing
// a couple of functors (Has_on_negative_side_2, for example), which are necessary in
// CGAL::Min_quadrilateral_default_traits_2. And while we know we have a generic case here,
// it's a bit tedious to get something well-defined in the generic case: for example,
// what should Has_on_negative_side_2 do if the Line_3 is orthogonal to the projection plane?
//
// So, easier to just bail out to real 2D...
template <typename Traits, typename PointRange>
std::pair<typename Traits::FT, typename Traits::FT>
compute_2D_deviation(const PointRange& points,
                     const PROJECTION_DIRECTION dir,
                     const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Traits::Point_2                                   Point_2;
  typedef typename Traits::Vector_2                                  Vector_2;
  typedef typename Traits::Point_3                                   Point_3;

  std::vector<Point_2> points_2D;
  points_2D.reserve(points.size());
  for(const Point_3& pt : points)
  {
    if(dir == ALONG_X)
      points_2D.emplace_back(pt.y(), pt.z());
    else if(dir == ALONG_Y)
      points_2D.emplace_back(pt.x(), pt.z());
    else if(dir == ALONG_Z)
      points_2D.emplace_back(pt.x(), pt.y());
  }

  std::vector<Point_2> extreme_points;
  ch_akl_toussaint(points_2D.begin(), points_2D.end(), std::back_inserter(extreme_points), traits);

  CGAL::Polygon_2<Traits> pol;
  CGAL::Min_quadrilateral_default_traits_2<Traits> mrt;
  CGAL::min_rectangle_2(extreme_points.begin(), extreme_points.end(), std::back_inserter(pol), mrt);

  if(pol.size() == 4 || !pol.is_simple() || pol.is_clockwise_oriented())
    return std::make_pair(0., 0.);

  // Compute the angle between the angle necessary to rotate the rectangle onto the reference frame
  auto bot_pos = pol.bottom_vertex();
  auto next_pos = bot_pos;
  ++next_pos;
  if(next_pos == pol.vertices_end())
    next_pos = pol.begin();

  const Point_2& p = *bot_pos;
  const Point_2& q = *next_pos;

  const Vector_2 pq = traits.construct_vector_2_object()(p, q);
  double n = sqrt(to_double(traits.compute_squared_length_2_object()(pq)));

  if(n == 0.) // degenerate input, maybe? Let's just not do anything
    return std::make_pair(pol.area(), 0.);

  const double dot = pq.x(); // that's the scalar product of PQ with V(1, 0) (Ox)
  double cosine = dot / n;

  if(cosine > 1.)
    cosine = 1.;
  if(cosine < -1.)
    cosine = -1.;

  double theta = std::acos(cosine);
  if(theta > 0.25 * CGAL_PI) // @todo is there a point to this
    theta = 0.5 * CGAL_PI - theta;

  return std::make_pair(pol.area(), FT{theta});
}

template <typename PointRange, typename Traits>
void optimize_along_OBB_axes(typename Traits::Matrix& rot,
                             const PointRange& points,
                             const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Traits::Point_3                                   Point;
  typedef typename Traits::Matrix                                    Matrix;
  typedef typename Traits::Vector                                    Vector;

  CGAL_static_assertion((std::is_same<typename boost::range_value<PointRange>::type, Point>::value));

  std::vector<Point> rotated_points;
  rotated_points.reserve(points.size());

  FT xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = FT{(std::numeric_limits<double>::max)()};
  xmax = ymax = zmax = FT{std::numeric_limits<double>::lowest()};

  for(const Point& pt : points)
  {
    Vector pv(3);
    pv.set(0, pt.x());
    pv.set(1, pt.y());
    pv.set(2, pt.z());
    pv = rot * pv;

    rotated_points.emplace_back(pv(0), pv(1), pv(2));

    xmin = (std::min)(xmin, pv(0));
    ymin = (std::min)(ymin, pv(1));
    zmin = (std::min)(zmin, pv(2));
    xmax = (std::max)(xmax, pv(0));
    ymax = (std::max)(ymax, pv(1));
    zmax = (std::max)(zmax, pv(2));
  }

  const FT lx = xmax - xmin;
  const FT ly = ymax - ymin;
  const FT lz = zmax - zmin;

  std::array<FT, 3> angles;
  std::array<FT, 3> volumes;

  FT area_xy;
  std::tie(area_xy, angles[0]) = compute_2D_deviation(rotated_points, ALONG_Z, traits);
  volumes[0] = lz * area_xy;

  FT area_xz;
  std::tie(area_xz, angles[1]) = compute_2D_deviation(rotated_points, ALONG_Y, traits);
  volumes[1] = ly * area_xz;

  FT area_yz;
  std::tie(area_yz, angles[2]) = compute_2D_deviation(rotated_points, ALONG_X, traits);
  volumes[2] = lx * area_yz;

  auto it = std::min_element(volumes.begin(), volumes.end());
  typename std::iterator_traits<decltype(it)>::difference_type d = std::distance(volumes.begin(), it);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
  std::cout << "volumes: " << volumes[0] << " " << volumes[1] << " " << volumes[2] << std::endl;
  std::cout << "angles: " << angles[0] << " " << angles[1] << " " << angles[2] << std::endl;
  std::cout << "min at " << d << std::endl;
#endif

  if(d == 0) // Along_Z
  {
    const double c = std::cos(angles[0]);
    const double s = std::sin(angles[0]);

    Matrix opt;
    opt.set(0, 0,  c); opt.set(0, 1, s); opt.set(0, 2, 0);
    opt.set(1, 0, -s); opt.set(1, 1, c); opt.set(1, 2, 0);
    opt.set(2, 0,  0); opt.set(2, 1, 0); opt.set(2, 2, 1);

    rot = opt * rot;
  }
  else if(d == 1) // Along_Y
  {
    const double c = std::cos(angles[1]);
    const double s = std::sin(angles[1]);

    Matrix opt;
    opt.set(0, 0, c); opt.set(0, 1, 0); opt.set(0, 2, -s);
    opt.set(1, 0, 0); opt.set(1, 1, 1); opt.set(1, 2, 0);
    opt.set(2, 0, s); opt.set(2, 1, 0); opt.set(2, 2, c);

    rot = opt * rot;
  }
  else if(d == 2) // Along_X
  {
    const double c = std::cos(angles[2]);
    const double s = std::sin(angles[2]);

    Matrix opt;
    opt.set(0, 0, 1); opt.set(0, 1,  0); opt.set(0, 2, 0);
    opt.set(1, 0, 0); opt.set(1, 1,  c); opt.set(1, 2, s);
    opt.set(2, 0, 0); opt.set(2, 1, -s); opt.set(2, 2, c);

    rot = opt * rot;
  }
  else
  {
    CGAL_assertion(false);
  }
}

// This operation makes no sense if an exact number type is used, so skip it, if so
template <typename Traits,
          typename IsFTExact = typename Algebraic_structure_traits<typename Traits::FT>::Is_exact>
struct Optimizer_along_axes
{
  template <typename PointRange>
  void operator()(typename Traits::Matrix& rot, const PointRange& points, const Traits& traits)
  {
    return optimize_along_OBB_axes(rot, points, traits);
  }
};

template <typename Traits>
struct Optimizer_along_axes<Traits, CGAL::Tag_true>
{
  template <typename PointRange>
  void operator()(typename Traits::Matrix&, const PointRange&, const Traits&) { }
};

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_OPTIMIZE_2_H
