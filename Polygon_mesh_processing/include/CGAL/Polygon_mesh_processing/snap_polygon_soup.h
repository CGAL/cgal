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

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_ROUNDING_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_ROUNDING_H

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

namespace CGAL
{

namespace Polygon_mesh_processing
{

/**
* 
*
* rounds the coordinates of the points so that they fit in doubles while keeping the model intersection free.
* The input can be any triangle soup and the output is an intersection-free triangle soup with Hausdorff distance
* between the input and the output bounded by M*2^-gs*k where M is the maximum absolute coordinate in the model, gs the snap_grid_size (see below) and k the number of iteration
* performed by the algorithm.
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam TriangleRange a model of the concepts `RandomAccessContainer`, `BackInsertionSequence` and `Swappable`, whose
* value type is a model of the concept `RandomAccessContainer` whose value type is convertible to `std::size_t` and that
* is constructible from an `std::initializer_list<std::size_t>` of size 3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param soup_triangles each element in the range describes a triangle using the indexed position of the points in `soup_points`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
*     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
*     \cgalParamDefault{`CGAL::Sequential_tag`}
*   \cgalParamNEnd
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `soup_points`}
*     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
*   \cgalParamNEnd
*   \cgalParamNBegin{snap_grid_size}
*     \cgalParamDescription{Scale the points to [-2^gs, 2^gs] where gs is the snap_grid_size before to round them on integer.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{23}
*     \cgalParamExtra{Must be lower than 52.}
*   \cgalParamNEnd
*   \cgalParamNBegin{numbers_of_iteration}
*     \cgalParamDescription{Maximum number of iteration performed by the algorithm.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{20}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \return `true` if the modified triangle soup is free from self-intersection, and `false` if the algorithm was not 
* able to provide such a triangle soup within the number of iterations.
*/
template <typename PointRange, typename PolygonRange, class NamedParameters = parameters::Default_named_parameters>
bool snap_polygon_soup(PointRange &points,
                       PolygonRange &triangles,
                       const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetPolygonSoupGeomTraits<PointRange, NamedParameters>::type GT;
  typedef typename GetPointMap<PointRange, NamedParameters>::const_type    Point_map;
  Point_map pm = choose_parameter<Point_map>(get_parameter(np, internal_np::point_map));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

#ifndef CGAL_LINKED_WITH_TBB
  static_assert (!parallel_execution,
                 "Parallel_tag is enabled but TBB is unavailable.");
#endif

  using Point_3 = std::remove_cv_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Kernel = typename Kernel_traits<Point_3>::Kernel;

  //Get the grid size from the named parameter, the grid size could not be greater than 52
  const unsigned int grid_size = (std::min)(52,choose_parameter(get_parameter(np, internal_np::snap_grid_size), 23));
  const unsigned int max_nb_of_iteration = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 20);

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
  std::cout << "Compute the scaling of the coordinates" << std::endl;
#endif

  //
  auto exp = [](const double v)
  {
    int n;
    frexp(v, &n);
    return n;
  };

  Bbox_3 bb = bbox_3(points.begin(), points.end());
  std::array<double, 3> max_abs{(std::max)(-bb.xmin(), bb.xmax()),
                                (std::max)(-bb.ymin(), bb.ymax()),
                                (std::max)(-bb.zmin(), bb.zmax())};
  std::array<double, 3> scale{std::pow(2, grid_size - exp(max_abs[0]) - 1),
                              std::pow(2, grid_size - exp(max_abs[1]) - 1),
                              std::pow(2, grid_size - exp(max_abs[2]) - 1)};

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

  for (std::size_t i = 0; i <= max_nb_of_iteration; ++i)
  {
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Start Iteration " << i << std::endl;
    std::cout << "Model size: " << points.size() << " " << triangles.size() << std::endl;
#endif

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Round all coordinates on doubles" << std::endl;
#endif
    for (Point_3 &p : points)
      p = Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));
    repair_polygon_soup(points, triangles, np);

    std::vector<std::pair<std::size_t, std::size_t>> pairs_of_intersecting_triangles;
    triangle_soup_self_intersections(points, triangles, std::back_inserter(pairs_of_intersecting_triangles));

    if (pairs_of_intersecting_triangles.empty())
    {
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "End of the snapping" << std::endl;
#endif
      CGAL_assertion(!does_triangle_soup_self_intersect(points, triangles));
      return true;
    }

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Number of pairs of intersecting triangles: " << pairs_of_intersecting_triangles.size() << std::endl;
#endif

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Snap the coordinates of the vertices of the intersecting triangles" << std::endl;
#endif
    //Get all the snap version of the points of the vertices of the intersecting triangles
    //Note: points will not be modified here, they will be modified in the next for loop

#if 1
    std::vector<Point_3> snap_points;
    snap_points.reserve(pairs_of_intersecting_triangles.size() * 3);

    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
        snap_points.emplace_back(snap_p(points[pi]));
      for (size_t pi : triangles[pair.second])
        snap_points.emplace_back(snap_p(points[pi]));
    }

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Snap the coordinates of the vertices close-by the previous ones" << std::endl;
#endif
    //TODO: TBB version of this for loop

    std::sort(snap_points.begin(), snap_points.end());
    snap_points.erase(std::unique(snap_points.begin(), snap_points.end()), snap_points.end());

    //If the snapped version of a point correspond to one of the previous point, we snap it too
    for (Point_3 &p : points)
    {
      Point_3 p_snap = snap_p(p);
      if (std::binary_search(snap_points.begin(), snap_points.end(), p_snap))
        p = p_snap;
    }
#else
    std::map<Point_3, size_t> snap_points;
    std::size_t index=0;
    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
      for (size_t pi : triangles[pair.second])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
    }

    //TODO: TBB version of this for loop

    std::vector<std::vector<size_t>> identical_points(index);
    //If the snapped version of a point correspond to one of the previous point, we snap it too
    for (size_t i=0; i!=points.size(); ++i)
    {
      Point_3 p_snap = snap_p(points[i]);
      auto it=snap_points.find(p_snap);
      if (it!=snap_points.end()){
        identical_points[it->second].push_back(i);
      }
    }

    for(const auto &v: identical_points){
      if(v.size()>1){
        std::array<double, 3> a{0,0,0};
        for(auto i: v){
          a[0]+=points[i].x();
          a[1]+=points[i].y();
          a[2]+=points[i].z();
        }
        a[0]/=v.size();
        a[1]/=v.size();
        a[2]/=v.size();
        for(auto i: v){
          points[i]=Point_3(a[0],a[1],a[2]);
        }
      }
    }
#endif


    repair_polygon_soup(points, triangles, np);
    //TODO do not pass all triangles
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Autorefine the soup" << std::endl;
    std::cout << "Model size: " << points.size() << " " << triangles.size() << std::endl;
#endif
    autorefine_triangle_soup(points, triangles, np);
  }
  return false;
}

} } //end of CGAL::Polygon_mesh_processing namespace

#endif //CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_ROUNDING_H
