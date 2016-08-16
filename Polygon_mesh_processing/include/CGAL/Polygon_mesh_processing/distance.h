// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
#define CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H

#include <algorithm>
#include <cmath>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/utility.h>
#include <boost/foreach.hpp>
//#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
//#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/spatial_sort.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <CGAL/atomic.h>
#endif // CGAL_LINKED_WITH_TBB


namespace CGAL{
namespace Polygon_mesh_processing {
namespace PMP = CGAL::Polygon_mesh_processing;
namespace internal{
template <class Kernel, class OutputIterator>
OutputIterator
triangle_grid_sampling( const typename Kernel::Point_3& p0,
                        const typename Kernel::Point_3& p1,
                        const typename Kernel::Point_3& p2,
                        double distance,
                        OutputIterator out)
{
  const double d_p0p1 = CGAL::sqrt( CGAL::squared_distance(p0, p1) );
  const double d_p0p2 = CGAL::sqrt( CGAL::squared_distance(p0, p2) );

  const double n = (std::max)(std::ceil( d_p0p1 / distance ),
                              std::ceil( d_p0p2 / distance ));

  for (double i=1; i<n; ++i)
    for (double j=1; j<n-i; ++j)
    {
      const double c0=(1-(i+j)/n), c1=i/n, c2=j/n;
      *out++=typename Kernel::Point_3(
              p0.x()*c0+p1.x()*c1+p2.x()*c2,
              p0.y()*c0+p1.y()*c1+p2.y()*c2,
              p0.z()*c0+p1.z()*c1+p2.z()*c2
            );
    }
  return out;
}

template <class Kernel, class OutputIterator>
OutputIterator
triangle_grid_sampling(const typename Kernel::Triangle_3& t, double distance, OutputIterator out)
{
  return triangle_grid_sampling<Kernel>(t[0], t[1], t[2], distance, out);
}

} //end of namespace internal

template <class Kernel, class TriangleRange, class OutputIterator>
OutputIterator
sample_triangles(const TriangleRange& triangles, double distance, OutputIterator out)
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  std::set< std::pair<Point_3, Point_3> > sampled_edges;

  // sample edges but skip endpoints
  BOOST_FOREACH(const Triangle_3& t, triangles)
  {
    for (int i=0;i<3; ++i)
    {
      const Point_3& p0=t[i];
      const Point_3& p1=t[(i+1)%3];
      if ( sampled_edges.insert(CGAL::make_sorted_pair(p0, p1)).second )
      {
        const double d_p0p1 = CGAL::sqrt( CGAL::squared_distance(p0, p1) );

        const double nb_pts = std::ceil( d_p0p1 / distance );
        const Vector_3 step_vec = (p1 - p0) / nb_pts;
        for (double i=1; i<nb_pts; ++i)
        {
          *out++=p0 + step_vec * i;
        }
      }
    }
  }

  // sample triangles
  BOOST_FOREACH(const Triangle_3& t, triangles)
    out=internal::triangle_grid_sampling<Kernel>(t, distance, out);

  //add endpoints
  std::set< Point_3 > endpoints;
  BOOST_FOREACH(const Triangle_3& t, triangles)
    for(int i=0; i<3; ++i)
    {
      if ( endpoints.insert(t[i]).second ) *out++=t[i];
    }
  return out;
}

#ifdef CGAL_LINKED_WITH_TBB
template <class AABB_tree, class Point_3>
struct Distance_computation{
  const AABB_tree& tree;
  const std::vector<Point_3>& sample_points;
  cpp11::atomic<double>* distance;

  Distance_computation(const AABB_tree& tree, const std::vector<Point_3>& sample_points, cpp11::atomic<double>* d)
    : tree(tree)
    , sample_points(sample_points)
    , distance(d)
  {}

  void
  operator()(const tbb::blocked_range<std::size_t>& range) const
  {
    double hdist = 0;
    Point_3 hint = sample_points.front();

    for( std::size_t i = range.begin(); i != range.end(); ++i)
    {
      hint = tree.closest_point(sample_points[i], hint);
      double d = CGAL::sqrt( squared_distance(hint,sample_points[i]) );
      if (d>hdist) hdist=d;
    }

    if (hdist > distance->load(cpp11::memory_order_acquire))
      distance->store(hdist, cpp11::memory_order_release);
  }
};
#endif

template <class Concurrency_tag, class Kernel, class TriangleMesh, class VertexPointMap = typename boost::property_map<TriangleMesh,
                                                                                                                       CGAL::vertex_point_t>::type>
double approximated_Hausdorff_distance(
  std::vector<typename Kernel::Point_3>& sample_points,
  TriangleMesh& m,
  std::size_t nb_sample_points)
{
 typedef Point_3<Kernel> Point_3;
 bool is_triangle = is_triangle_mesh(m);
  CGAL_assertion_msg (is_triangle,
        "Mesh is not triangulated. Distance computing impossible.");
  Random_points_in_triangle_mesh_3<TriangleMesh, VertexPointMap>
      g(m);
  CGAL::cpp11::copy_n(g, nb_sample_points, std::back_inserter(sample_points));
  #ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "Nb sample points " << sample_points.size() << "\n";
  #endif
  spatial_sort(sample_points.begin(), sample_points.end());

  typedef typename TriangleMesh::face_iterator Triangle_iterator;
  typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef AABB_traits<Kernel, Primitive> Traits;
  typedef AABB_tree< Traits > Tree;

  Tree tree( faces(m).first, faces(m).second, m);
  tree.accelerate_distance_queries();
  tree.build();

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<Concurrency_tag,Parallel_tag>::value)
  {
    cpp11::atomic<double> distance(0);
    Distance_computation<Tree, typename Kernel::Point_3> f(tree, sample_points, &distance);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, sample_points.size()), f);
    return distance;
  }
  else
#endif
  {
    double hdist = 0;
    typename Traits::Point_3 hint = sample_points.front();
    BOOST_FOREACH(const typename Traits::Point_3& pt, sample_points)
    {
      hint = tree.closest_point(pt, hint);
      typename Kernel::FT dist = squared_distance(hint,pt);
      double d = to_double(CGAL::approximate_sqrt(dist));
      if (d>hdist) hdist=d;
    }
    return hdist;
  }
}

template <class Concurrency_tag, class Kernel, class TriangleMesh,
          class VertexPointMap1 = typename boost::property_map<TriangleMesh,
                                                               CGAL::vertex_point_t>::type,
          class VertexPointMap2 = typename boost::property_map<TriangleMesh,
                                                               CGAL::vertex_point_t>::type>
double approximated_Hausdorff_distance(
   TriangleMesh& m1,
   TriangleMesh& m2,
  int  nb_points
)
{
  std::vector<typename Kernel::Point_3> sample_points;
  Random_points_in_triangle_mesh_3<TriangleMesh, VertexPointMap1>
      g(m1);
  CGAL::cpp11::copy_n(g, nb_points, std::back_inserter(sample_points));
  return approximated_Hausdorff_distance<Concurrency_tag, Kernel, TriangleMesh, VertexPointMap2>(sample_points, m2,4000);
}

template <class Concurrency_tag, class Kernel, class TriangleMesh, class VertexPointMap1 = typename boost::property_map<TriangleMesh,
                                                                                                                       CGAL::vertex_point_t>::type,
          class VertexPointMap2 = typename boost::property_map<TriangleMesh,
                                                                                            CGAL::vertex_point_t>::type>
double approximated_symmetric_Hausdorff_distance(
  TriangleMesh& m1,
  TriangleMesh& m2,
 int  nb_points
)
{
  return (std::max)(
    approximated_Hausdorff_distance<Concurrency_tag, Kernel, TriangleMesh, VertexPointMap1, VertexPointMap2>(m1, m2, nb_points),
    approximated_Hausdorff_distance<Concurrency_tag, Kernel, TriangleMesh, VertexPointMap2, VertexPointMap1>(m2, m1, nb_points)
  );
}

/// \todo test different strategies and put the better one in `approximated_Hausdorff_distance()`
/// for particular cases one can still use a specific sampling method together with `max_distance_to_triangle_mesh()`

enum Sampling_method{
 RANDOM_UNIFORM =0, /**< points are generated in a random and uniform way, depending on the area of each triangle.*/
 GRID,/**< points are generated in a grid, with a minimum of one point per triangle.*/
 MONTE_CARLO /**< points are generated randomly in each triangle, proportionally to the face area with a minimum
               * of 1 pt per triangle.*/
};
/** fills `sampled_points` with points taken on the mesh in a manner depending on `method`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph` that has an internal property map
 *         for `CGAL::vertex_point_t`
 * @param m the triangle mesh that will be sampled
 * @param precision depends on the value of `method` :
 * case RANDOM_UNIFORM : the number of points per squared area unit
 * case GRID : the distance between the points
 * case MONTE_CARLO : the number of points per squared area unit
 * @tparam method a Sampling_method.
 * @tparam Sampling_method defines the method of sampling.
 *                         Possible values are `RANDOM_UNIFORM`,
 *                         and `GRID` and `MONTE_CARLO.
 */
template<class Kernel, class TriangleMesh>
void sample_triangle_mesh(const TriangleMesh& m,
                          double precision,
                          std::vector<typename Kernel::Point_3>& sampled_points,
                          Sampling_method method = RANDOM_UNIFORM)
{
 switch(method)
 {
 case RANDOM_UNIFORM:
 {
  std::size_t nb_points =  std::ceil(precision * PMP::area(m,
                                                           PMP::parameters::geom_traits(Kernel())));
   Random_points_in_triangle_mesh_3<TriangleMesh>
      g(m);
   CGAL::cpp11::copy_n(g, nb_points, std::back_inserter(sampled_points));
   return;
 }
 case GRID:
 {
   typedef typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type Pmap;
  Pmap pmap = get(vertex_point, m);
   BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(m))
   {
     typename Kernel::Point_3 points[3];
     typename TriangleMesh::Halfedge_around_face_circulator hc(halfedge(f,m), m);
     for(int i=0; i<3; ++i)
     {
       points[i] = get(pmap, target(*hc, m));
       ++hc;
    }
    internal::triangle_grid_sampling<Kernel>(points[0], points[1], points[2],precision, std::back_inserter(sampled_points));
  }
  return;
 }
 case MONTE_CARLO:
  std::size_t nb_points =  std::ceil(precision * PMP::area(m,
                                                           PMP::parameters::geom_traits(Kernel())));
  typedef typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type Pmap;
 Pmap pmap = get(vertex_point, m);
  std::vector<typename Kernel::Triangle_3> triangles;
  BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(m))
  {
  //create the triangles and store them
    typename Kernel::Point_3 points[3];
    typename TriangleMesh::Halfedge_around_face_circulator hc(halfedge(f,m), m);
    for(int i=0; i<3; ++i)
    {
      points[i] = get(pmap, target(*hc, m));
      ++hc;
    }
    triangles.push_back(typename Kernel::Triangle_3(points[0], points[1], points[2]));
    //sample a single point in all triangles(to have at least 1 pt/triangle)
    Random_points_in_triangle_3<typename Kernel::Point_3> g(points[0], points[1], points[2]);
    CGAL::cpp11::copy_n(g, 1, std::back_inserter(sampled_points));
  }
  //sample the triangle range uniformly
  Random_points_in_triangles_3<typename Kernel::Point_3 >
    g(triangles);
 CGAL::cpp11::copy_n(g, nb_points, std::back_inserter(sampled_points));
  return;
 }
}
/// \todo add a plugin in the demo to display the distance between 2 meshes as a texture (if not complicated)

template< class Concurrency_tag,
          class Kernel,
          class TriangleMesh,
          class PMap1,
          class PMap2>
double approximated_Hausdorff_distance( TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        double precision,
                                        const PMap1&,
                                        const PMap2&)
{
 std::size_t nb_points =  std::max(std::ceil(to_double(precision * PMP::area(tm1,
                                                                 PMP::parameters::geom_traits(Kernel())))),
                                   std::ceil(to_double(precision * PMP::area(tm2,
                                                                 PMP::parameters::geom_traits(Kernel())))));
 return approximated_Hausdorff_distance<Concurrency_tag,Kernel,TriangleMesh, PMap1, PMap2>(tm1, tm2, nb_points);
}
// documented functions
/**
 * \ingroup PMP_distance_grp
 * computes the approximated Hausdorff distance of `tm1` from `tm2` by
 * generating a uniform random point sampling on `tm1`, and by then
 * returning the distance of the furthest point from `tm2`.
 *
 * A parallel version is provided and requires the executable to be
 * linked against the <a href="http://www.threadingbuildingblocks.org">Intel TBB library</a>.
 * To control the number of threads used, the user may use the `tbb::task_scheduler_init` class.
 * See the <a href="http://www.threadingbuildingblocks.org/documentation">TBB documentation</a>
 * for more details.
 *
 * @tparam Concurrency_tag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag`
 *                         and `Parallel_tag`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph` that has an internal property map
 *         for `CGAL::vertex_point_t`
 * @tparam NamedParameters1 a sequence of \ref namedparameters for `tm1`
 * @tparam NamedParameters2 a sequence of \ref namedparameters for `tm2`
 *
 * @param tm1 the triangle mesh that will be sampled
 * @param tm2 the triangle mesh to compute the distance to
 * @param precision the number of points per squared area unit
 * @param np1 optional sequence of \ref namedparameters for `tm1` among the ones listed below
 * @param np2 optional sequence of \ref namedparameters for `tm2` among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 * \cgalNamedParamsEnd
 * \todo When using TBB, runtime is slower. The chunk size should probably be tuned. see PMP/test/test_pmp_distance.cpp
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximated_Hausdorff_distance( TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        double precision,
                                        const NamedParameters1& np1,
                                        const NamedParameters2& np2)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters1>::type Geom_traits;

  return approximated_Hausdorff_distance<Concurrency_tag, Geom_traits>(
    tm1, tm2,
    precision,
    choose_const_pmap(get_param(np1, boost::vertex_point),
                                tm1,
                                vertex_point),
    choose_const_pmap(get_param(np2, boost::vertex_point),
                                tm2,
                                vertex_point)

  );
}

/**
 * \ingroup PMP_distance_grp
 * computes the approximated symmetric Hausdorff distance between `tm1` and `tm2`.
 * It returns the maximum of `approximated_Hausdorff_distance(tm1,tm2)` and
 * `approximated_Hausdorff_distance(tm1,tm2)`.
 *
 * \copydetails CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance()
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximated_symmetric_Hausdorff_distance(
  TriangleMesh& tm1,
  TriangleMesh& tm2,
  double precision,
  const NamedParameters1& np1,
  const NamedParameters2& np2)
{
  return (std::max)(
    approximated_Hausdorff_distance<Concurrency_tag>(tm1,tm2,precision,np1,np2),
    approximated_Hausdorff_distance<Concurrency_tag>(tm2,tm1,precision,np2,np1)
  );
}

/// \todo document and implement me
/// \todo find a way to define precision through named parameters
/**
 * \ingroup PMP_distance_grp
 * computes the approximated Hausdorff distance between `points` and `tm`.
 * \copydetails CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance()
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class PointRange,
          class NamedParameters>
double max_distance_to_triangle_mesh(const PointRange& points,
                                     TriangleMesh& tm,
                                     double precision,
                                     const NamedParameters& np)
{
 typedef typename GetGeomTraits<TriangleMesh,
                                NamedParameters>::type Geom_traits;
 std::vector<typename PointRange::value_type> sample_points;
 BOOST_FOREACH(typename PointRange::value_type point, points)
   sample_points.push_back(point);

 std::size_t nb_points =  std::ceil(to_double(precision * PMP::area(tm,
                                                          PMP::parameters::geom_traits(Geom_traits()))));
  return approximated_Hausdorff_distance<Concurrency_tag, Geom_traits, TriangleMesh/*,
     choose_const_pmap(get_param(np, boost::vertex_point),
                       tm,
                       vertex_point)*/>
     (sample_points,tm, nb_points);
}

/// \todo document and implement me
///       see with @sgiraudot for the implementation
///       that should be put in a seperate header file
///       with copyright INRIA
///       Maybe find a better name too
template< class Concurrency_tag,
          class TriangleMesh,
          class PointRange,
          class NamedParameters>
double max_distance_to_point_set(TriangleMesh& tm,
                                 const PointRange& points,
                                 const NamedParameters& np)
{
  return 0;
}

// convenience functions with default parameters

template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1>
double approximated_Hausdorff_distance( TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        double precision,
                                        const NamedParameters1& np1)
{
  return approximated_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision, np1, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double approximated_Hausdorff_distance( TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        double precision)
{
  return approximated_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision,
    parameters::all_default(),
    parameters::all_default());
}


template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1>
double approximated_symmetric_Hausdorff_distance(TriangleMesh& tm1,
                                                 TriangleMesh& tm2,
                                                 double precision,
                                                 const NamedParameters1& np1)
{
  return approximated_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision, np1, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1>
double approximated_symmetric_Hausdorff_distance(TriangleMesh& tm1,
                                                 TriangleMesh& tm2,
                                                 double precision)
{
  return approximated_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision,
    parameters::all_default(),
    parameters::all_default());
}

} } // end of namespace CGAL::Polygon_mesh_processing


#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
