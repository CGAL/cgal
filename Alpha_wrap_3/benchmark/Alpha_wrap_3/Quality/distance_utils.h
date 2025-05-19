// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Michael Hemmer
//                 Cedric Portaneri

#ifndef CGAL_ALPHA_WRAP_3_BENCHMARK_ALPHA_WRAP_3_QUALITY_DISTANCE_H_
#define CGAL_ALPHA_WRAP_3_BENCHMARK_ALPHA_WRAP_3_QUALITY_DISTANCE_H_

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Search_traits_3.h>

namespace Aw3i {

enum Distance_metric { HAUSDORFF = 0, MEAN = 1, RMS = 2 };

template <typename Point, typename AABBTree>
inline double approximate_hausdorff_distance(const std::vector<Point>& sample_points,
                                             const AABBTree& tree,
                                             Point& hint)
{
  double hdist = 0;
  for(const Point& pt : sample_points)
  {
    hint = tree.closest_point(pt, hint);
    auto dist = CGAL::squared_distance(hint, pt);
    double d = CGAL::to_double(CGAL::approximate_sqrt(dist));
    if(d > hdist)
      hdist = d;
  }

  return hdist;
}

template <typename Point, typename AABBTree>
inline double approximate_mean_distance(const std::vector<Point>& sample_points,
                                        const AABBTree& tree,
                                        Point& hint)
{
  double mdist = 0;
  for(const Point& pt : sample_points)
  {
    hint = tree.closest_point(pt, hint);
    auto dist = CGAL::squared_distance(hint, pt);
    double d = CGAL::to_double(CGAL::approximate_sqrt(dist));
    mdist += d;
  }

  return mdist / sample_points.size();
}

template <typename Point, typename AABBTree>
inline double approximate_rms_distance(const std::vector<Point>& sample_points,
                                       const AABBTree& tree,
                                       Point& hint)
{
  double rmsdist = 0;
  for(const Point& pt : sample_points)
  {
    hint = tree.closest_point(pt, hint);
    auto dist = CGAL::squared_distance(hint, pt);
    rmsdist += CGAL::to_double(dist);
  }

  return CGAL::to_double(CGAL::approximate_sqrt(rmsdist / sample_points.size()));
}

template <typename TriangleMesh>
inline double approximate_distance(const TriangleMesh& tm1,
                                   const TriangleMesh& tm2,
                                   const Distance_metric& metric)
{
  using GT = typename CGAL::GetGeomTraits<TriangleMesh>::type;
  using Point_3 = typename GT::Point_3;

  using Primitive = CGAL::AABB_face_graph_triangle_primitive<TriangleMesh>;
  using AABB_traits = CGAL::AABB_traits_3<GT, Primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  std::vector<Point_3> original_sample_points;
  CGAL::Polygon_mesh_processing::sample_triangle_mesh(tm1, std::back_inserter(original_sample_points),
                                                      CGAL::parameters::all_default());

  std::vector<Point_3> sample_points(std::begin(original_sample_points),
                                     std::end(original_sample_points));
  CGAL::spatial_sort(sample_points.begin(), sample_points.end());

  AABB_tree tree(faces(tm2).first, faces(tm2).second, tm2);
  tree.build();

  auto vpm_2 = get(CGAL::vertex_point, tm2);
  Point_3 hint = get(vpm_2, *vertices(tm2).first);

  if(metric == HAUSDORFF)
    return approximate_hausdorff_distance(sample_points, tree, hint);
  else if(metric == MEAN)
    return approximate_mean_distance(sample_points, tree, hint);
  else if(metric == RMS)
    return approximate_rms_distance(sample_points, tree, hint);
  else
    std::cerr << "Metric unknown\n" << std::endl;

  return -1.0;
}

template <typename TriangleMesh>
double get_longest_diag_bbox(const TriangleMesh& tm)
{
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(tm);
  return std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                   CGAL::square(bbox.ymax() - bbox.ymin()) +
                   CGAL::square(bbox.zmax() - bbox.zmin()));
}

template <typename TriangleMesh>
inline double approximate_distance_relative_to_bbox(const TriangleMesh& tm1,
                                                    const TriangleMesh& tm2,
                                                    const Distance_metric& metric)
{
  double longest_diag_length = get_longest_diag_bbox(tm1);
  return approximate_distance(tm1, tm2, metric) / longest_diag_length;
}

template <typename TriangleMesh, typename FT>
inline double approximate_distance_relative_to_bbox(const TriangleMesh& tm1,
                                                    const TriangleMesh& tm2,
                                                    const Distance_metric& metric,
                                                    const FT& longest_diag_length)
{
  return approximate_distance(tm1, tm2, metric) / CGAL::to_double(longest_diag_length);
}

} // namespace Aw3i

#endif // CGAL_CGAL_ALPHA_WRAP_3_BENCHMARK_ALPHA_WRAP_3_QUALITY_DISTANCE_H_
