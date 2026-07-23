// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//                 Bradley McCoy

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_POISSON_DISK_SAMPLING_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_POISSON_DISK_SAMPLING_H


#include <CGAL/point_generators_3.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#ifdef CGAL_USE_BSURF
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>
#endif

#include <vector>
#include <queue>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

enum Distance_version { EUCLIDEAN_DISTANCE, GEODESIC_DISTANCE };

template <class GeomTraits>
double euclideanDistancePoints(const typename GeomTraits::Point_3& source,
                               const typename GeomTraits::Point_3& target)
{
  return sqrt(CGAL::squared_distance(source,target));
}

#ifdef CGAL_USE_BSURF
template <class GeomTraits, class TriangleMesh>
double geodesiceApproximation(const Face_location<TriangleMesh, typename GeomTraits::FT>& source,
                              const Face_location<TriangleMesh, typename GeomTraits::FT>& target,
                              const TriangleMesh& mesh,
                              const Dual_geodesic_solver<double>& solver)
{
  std::vector<Edge_location<TriangleMesh, typename GeomTraits::FT>> edge_locations;
  CGAL::Polygon_mesh_processing::locally_shortest_path<double>(source, target, mesh, edge_locations, solver);

  return path_length<GeomTraits>(edge_locations,source,target,mesh);
}
#endif

//function to switch between geodesic and Euclidean distance
template <class GeomTraits, Distance_version V, class TriangleMesh>
double distancePoints(const TriangleMesh& /* mesh */,
                      const typename GeomTraits::Point_3& source,
                      const typename GeomTraits::Point_3& target,
                      const Face_location<TriangleMesh, typename GeomTraits::FT>& /* start */,
                      const Face_location<TriangleMesh, typename GeomTraits::FT>& /* end */
#ifdef CGAL_USE_BSURF
                      , const Dual_geodesic_solver<double>& solver
#endif
                      )
{
#ifdef CGAL_USE_BSURF
  if constexpr (V==GEODESIC_DISTANCE)
    return geodesiceApproximation<GeomTraits>(start, end, mesh, solver);
#endif
  if constexpr (V==EUCLIDEAN_DISTANCE)
    return euclideanDistancePoints<GeomTraits>(source, target);
  return 0;
}


template <class GeomTraits, Distance_version V, class TriangleMesh>
std::vector<typename boost::graph_traits<TriangleMesh>::face_descriptor>
faces_in_sub_mesh(const typename GeomTraits::Point_3& c,
                  typename boost::graph_traits<TriangleMesh>::face_descriptor fc,
                  const TriangleMesh& mesh,
                  double minDistance)
{
//  using Point = typename GeomTraits::Point_3;
  using Graph_traits = boost::graph_traits<TriangleMesh>;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;

  //Begin flooding
  face_descriptor fd = fc;

  std::vector<bool> selected(num_faces(mesh), false);
  std::vector<face_descriptor> selection;
  selected[fd] = true;
  selection.push_back(fd);

  auto do_queue_edge = [&](halfedge_descriptor h)
  {
    halfedge_descriptor hopp=opposite(h, mesh);
    if (is_border(hopp, mesh) || selected[face(hopp, mesh)]) return false;

    typename GeomTraits::Segment_3 edge(mesh.point(source(h,mesh)), mesh.point(target(h,mesh)));

//      return (distancePoints<V>(mesh, mesh.point(source(h,mesh)), c, tree)< 3*minDistance ||
//              distancePoints<V>(mesh, mesh.point(target(h,mesh)), c, tree)< 3*minDistance);


    return typename GeomTraits::Compare_squared_distance_3()(c, edge, sqrt(3*minDistance))!=CGAL::LARGER;
  };


  std::vector<halfedge_descriptor> queue;
  for (halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(fd, mesh), mesh))
    if (do_queue_edge(h))
      queue.push_back(opposite(h, mesh));

  while (!queue.empty())
  {
    halfedge_descriptor h = queue.back();
    face_descriptor f = face(h, mesh);
    queue.pop_back();
    if (!selected[f])
    {
      selected[f]=true;
      selection.push_back(f);
    }

    h=next(h, mesh);
    if (do_queue_edge(h)) queue.push_back(opposite(h, mesh));
    h=next(h, mesh);
    if (do_queue_edge(h)) queue.push_back(opposite(h, mesh));
  }

  return selection;
}

template <class GeomTraits, Distance_version V, class TriangleMesh, class PointsPerFaceMap>
bool
is_far_enough(const typename GeomTraits::Point_3 c,
              const Face_location<TriangleMesh, typename GeomTraits::FT> c_location,
              const TriangleMesh& mesh,
              double minDistance,
              std::vector<typename boost::graph_traits<TriangleMesh>::face_descriptor> selection,
              const PointsPerFaceMap& face_points
#ifdef CGAL_USE_BSURF
              , const Dual_geodesic_solver<double>& solver
#endif
              )
{
  for (typename boost::graph_traits<TriangleMesh>::face_descriptor f : selection)
  {
    for(const typename GeomTraits::Point_3& p  : face_points[f])
    {
      Face_location<TriangleMesh, typename GeomTraits::FT> p_location = locate_in_face(p, f, mesh);
      //Todo: Ask why  is distancePoints for Euclidean so much slower then just
      //calling euclideanDistancePoints?
      if (distancePoints<GeomTraits,V>(mesh, c, p, c_location, p_location
#ifdef CGAL_USE_BSURF
          , solver
#endif
                                      ) < minDistance)
      //if (euclideanDistancePoints(c, p) < minDistance)
      //if (geodesiceApproximation(c_location, p_location, mesh) < minDistance)
          return false;
    }
  }
  return true;
}

template <class GeomTraits, Distance_version V, class TriangleMesh>
std::vector<typename GeomTraits::Point_3>
poisson_disk_sampling(const TriangleMesh& mesh,
                      double minDistance,
                      std::size_t kMaxTries,
                      CGAL::Random& r)
{
  using Point = typename GeomTraits::Point_3;
  using Graph_traits = boost::graph_traits<TriangleMesh>;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using FT = typename GeomTraits::FT;
  using Face_location = Face_location<TriangleMesh, FT>;

#ifdef CGAL_USE_BSURF
  Dual_geodesic_solver<double> solver;
  if constexpr (V==GEODESIC_DISTANCE)
    init_geodesic_dual_solver(solver, mesh);
#endif

  std::vector<Point> points;
  std::queue<std::pair<Point, face_descriptor>> activePoints;

  Random_points_in_triangle_mesh_3<TriangleMesh> g(mesh, r);
  Point c = *g;
  face_descriptor fc = g.last_item_picked();

  std::vector<std::vector<Point>> face_points(num_faces(mesh));

  face_points[fc].push_back(c);
  activePoints.push(std::make_pair(c, fc));
  points.push_back(c);

  Point currentPoint;
  face_descriptor currentFace;
  while (!activePoints.empty())
  {
    std::tie(currentPoint, currentFace) = activePoints.front();
    activePoints.pop();

    std::vector<face_descriptor> selection = faces_in_sub_mesh<GeomTraits,V>(currentPoint, currentFace, mesh, minDistance);
    Face_location current_location = locate_in_face(currentPoint, currentFace, mesh);
    if (current_location.second[0]<0) current_location.second[0]=0;
    if (current_location.second[1]<0) current_location.second[1]=0;
    if (current_location.second[2]<0) current_location.second[2]=0;
    std::size_t k = 0;
    while (k < kMaxTries)
    {
      double angle = 2 * CGAL_PI * r.get_double();
      double distance = minDistance + minDistance * r.get_double();
      typename GeomTraits::Vector_2 dir(cos(angle),sin(angle));
      #warning find alternative
      CGAL_USE(distance);
      std::vector<Face_location> path /* = straightest_geodesic<GeomTraits>(current_location, dir, distance, mesh) */ ;
      Face_location newLocation = path.back();
      Point newPoint = construct_point(path.back(), mesh);
      if(is_far_enough<GeomTraits,V>(newPoint, newLocation, mesh, minDistance, selection, make_random_access_property_map(face_points)
#ifdef CGAL_USE_BSURF
        , solver
#endif
      ))
      {
        face_points[path.back().first].push_back(newPoint);
        activePoints.push(std::make_pair(newPoint,path.back().first));
        points.push_back(newPoint);
      }else
      {
        ++k;
      }
    }
  }

  return points;
}

} } } // CGAL::Polygon_mesh_processing::internal


#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_POISSON_DISK_SAMPLING_H
