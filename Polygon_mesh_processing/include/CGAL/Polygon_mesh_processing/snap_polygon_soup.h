// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Léo Valque

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_ROUNDING
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_ROUNDING

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

namespace CGAL
{

namespace Polygon_mesh_processing
{

template <typename PointRange, typename PolygonRange>
bool snap_polygon_soup(PointRange &points,
                                        PolygonRange &triangles)
{
  using Point_3 = std::remove_cv_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Kernel = typename Kernel_traits<Point_3>::Kernel;
  static constexpr int coarseness = 23;
  // TODO named parameters

  auto exp = [](const double v)
  {
    int n;
    frexp(v, &n);
    return n + 1;
  };

  Bbox_3 bb = bbox_3(points.begin(), points.end());
  std::array<double, 3> max_abs{(std::max)(-bb.xmin(), bb.xmax()),
                                (std::max)(-bb.ymin(), bb.ymax()),
                                (std::max)(-bb.zmin(), bb.zmax())};
  std::array<double, 3> scale{std::pow(2, coarseness - exp(max_abs[0])),
                              std::pow(2, coarseness - exp(max_abs[1])),
                              std::pow(2, coarseness - exp(max_abs[2]))};

  // If EPECK, use exact TODO
  auto snap = [](typename Kernel::FT x, double scale)
  {
    return std::ceil(CGAL::to_double(x * scale) + 0.5) / scale;
  };
  auto snap_p = [scale, snap](const Point_3 &p)
  {
    return Point_3(snap(p.x(),scale[0]),
                    snap(p.y(),scale[1]),
                    snap(p.z(),scale[2]) );
  };

  static constexpr std::size_t nb_max_of_iteration = 20; // TODO named parameter
  for (std::size_t i = 0; i <= nb_max_of_iteration; ++i)
  {
    std::cout << "after autorefine " << triangles.size() << std::endl;

    for (Point_3 &p : points)
      p = Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));
    repair_polygon_soup(points, triangles);

    std::vector<std::pair<std::size_t, std::size_t>> out;
    triangle_soup_self_intersections(points, triangles, std::back_inserter(out));

    if (out.empty())
    {
      CGAL_assertion(!does_triangle_soup_self_intersect(points, triangles));
      return true;
    }

    std::vector<Point_3> snap_points;
    snap_points.reserve(out.size() * 3);

    for (auto &pair : out)
    {
      for (size_t pi : triangles[pair.first])
        snap_points.emplace_back(snap_p(points[pi]));
      for (size_t pi : triangles[pair.second])
        snap_points.emplace_back(snap_p(points[pi]));
    }

    std::sort(snap_points.begin(), snap_points.end());
    snap_points.erase(std::unique(snap_points.begin(), snap_points.end()), snap_points.end());

    for (Point_3 &p : points)
    {
      Point_3 p_snap = snap_p(p);
      if (std::binary_search(snap_points.begin(), snap_points.end(), p_snap))
        p = p_snap;
    }

    repair_polygon_soup(points, triangles);
    //TODO do not pass all triangles
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Before autorefine " << triangles.size() << std::endl;
#endif
    autorefine_triangle_soup(points, triangles, parameters::concurrency_tag(Parallel_tag()));
  }
  return false;
}
}

}

#endif
