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
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/spatial_sort.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Polygon_mesh_processing/internal/mesh_to_point_set_hausdorff_distance.h>
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <CGAL/atomic.h>
#endif // CGAL_LINKED_WITH_TBB


namespace CGAL{
namespace Polygon_mesh_processing {
namespace internal{
template <class Kernel, class OutputIterator>
OutputIterator
triangle_grid_sampling( const typename Kernel::Point_3& p0,
                        const typename Kernel::Point_3& p1,
                        const typename Kernel::Point_3& p2,
                        double distance,
                        OutputIterator out)
{
    typename Kernel::Compute_squared_distance_3 squared_distance;
  const double d_p0p1 = to_double(approximate_sqrt( squared_distance(p0, p1) ));
  const double d_p0p2 = to_double(approximate_sqrt( squared_distance(p0, p2) ));

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

  // sample edges but      skip endpoints
  BOOST_FOREACH(const Triangle_3& t, triangles)
  {
    for (int i=0;i<3; ++i)
    {
      const Point_3& p0=t[i];
      const Point_3& p1=t[(i+1)%3];
      if ( sampled_edges.insert(CGAL::make_sorted_pair(p0, p1)).second )
      {
          typename Kernel::Compute_squared_distance_3 squared_distance;
        const double d_p0p1 = to_double(approximate_sqrt( squared_distance(p0, p1) ));

        const double nb_pts = std::ceil( d_p0p1 / distance );
        const Vector_3 step_vec =  typename Kernel::Construct_scaled_vector_3()(typename Kernel::Construct_vector_3()(p1, p0), typename Kernel::FT(1)/typename Kernel::FT(nb_pts));
        for (double i=1; i<nb_pts; ++i)
        {
          *out++=typename Kernel::Construct_translated_point_3()(p0, typename Kernel::Construct_scaled_vector_3()(step_vec , typename Kernel::FT(i)));
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
  Point_3 initial_hint;
  cpp11::atomic<double>* distance;

  Distance_computation(
          const AABB_tree& tree,
          const Point_3& p,
          const std::vector<Point_3>& sample_points,
          cpp11::atomic<double>* d)
    : tree(tree)
    , sample_points(sample_points)
    , initial_hint(p)
    , distance(d)
  {}

  void
  operator()(const tbb::blocked_range<std::size_t>& range) const
  {
    Point_3 hint = initial_hint;
    double hdist = 0;
    for( std::size_t i = range.begin(); i != range.end(); ++i)
    {
      hint = tree.closest_point(sample_points[i], hint);
      typename Kernel_traits<Point_3>::Kernel::Compute_squared_distance_3 squared_distance;
      double d = to_double(CGAL::approximate_sqrt( squared_distance(hint,sample_points[i]) ));
      if (d>hdist) hdist=d;
    }

    if (hdist > distance->load(cpp11::memory_order_acquire))
      distance->store(hdist, cpp11::memory_order_release);
  }
};
#endif

/**
 * \ingroup PMP_distance_grp
 * @brief enum used to select the sampling method in the function `sample_triangle_mesh()`
 */
enum Sampling_method{
  RANDOM_UNIFORM =0,
  GRID,
  MONTE_CARLO
};
/** \ingroup PMP_distance_grp
 * generates points taken on `tm` in a way depending on `method` and
 * outputs them to `out`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam OutputIterator a model of `OutputIterator` 
 *  holding objects of the same point type as
 *  the value type of the internal vertex point map of `tm`
 * @tparam Sampling_method defines the sampling method. Possible values are:
     - `RANDOM_UNIFORM`: points are generated in a random and uniform way, depending on the area of each triangle.
     - `GRID`: points are generated on a grid in each triangle, with a minimum of one point per triangle.
     - `MONTE_CARLO`: points are generated randomly in each triangle. The number of points per triangle is proportional to the triangle area with a minimum of 1.
 *
 * @param tm the triangle mesh that will be sampled
 * @param parameter depends on `method`:
     - `RANDOM_UNIFORM` and `MONTE_CARLO`: the number of points per squared area unit
     - `GRID`: the distance between two consecutive points in the grid
 *
 * @param out output iterator to be filled with sampled points
 * @param np an optional sequence of \ref namedparameters among the ones listed below
 * @param method defines the method of sampling
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map with the points associated to the vertices of `tm`. If this parameter is omitted,
 *    an internal property map for `CGAL::vertex_point_t` should be available for `TriangleMesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPDistanceTraits`\cgalParamEnd
 * \cgalNamedParamsEnd
 */
template<class OutputIterator, class TriangleMesh, class NamedParameters>
OutputIterator
sample_triangle_mesh(const TriangleMesh& tm,
                           double parameter,
                           OutputIterator out,
                           NamedParameters np,
                           Sampling_method method = RANDOM_UNIFORM)
{
    typedef typename GetGeomTraits<TriangleMesh,
            NamedParameters>::type Geom_traits;

    typedef typename GetVertexPointMap<TriangleMesh,
            NamedParameters>::const_type Vpm;

    Vpm pmap = choose_param(get_param(np, vertex_point),
                           get_const_property_map(vertex_point, tm));
    typedef Creator_uniform_3<typename Geom_traits::FT,
                              typename Geom_traits::Point_3> Creator;
  switch(method)
  {
    case RANDOM_UNIFORM:
    {
      std::size_t nb_points = std::ceil( to_double(
         area(tm, parameters::geom_traits(Geom_traits()))*parameter) );
      Random_points_in_triangle_mesh_3<TriangleMesh, Vpm, Creator>
          g(tm, pmap);
      CGAL::cpp11::copy_n(g, nb_points, out);
      return out;
    }
    case GRID:
    {
      std::vector<typename Geom_traits::Triangle_3> triangles;
      BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(tm))
      {
        //create the triangles and store them
        typename Geom_traits::Point_3 points[3];
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
        for(int i=0; i<3; ++i)
        {
          points[i] = get(pmap, target(hd, tm));
          hd = next(hd, tm);
        }
        triangles.push_back(typename Geom_traits::Triangle_3(points[0], points[1], points[2]));
      }
      sample_triangles<Geom_traits>(triangles, parameter, out);
      return out;
    }
    case MONTE_CARLO:
    {
      BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(tm))
      {
        std::size_t nb_points( (std::max)(
                    std::ceil(to_double(face_area(f,tm,parameters::geom_traits(Geom_traits())))*parameter),
                                         1.) );
        //create the triangles and store them
        typename Geom_traits::Point_3 points[3];
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
        for(int i=0; i<3; ++i)
        {
          points[i] = get(pmap, target(hd, tm));
          hd = next(hd, tm);
        }
        Random_points_in_triangle_3<typename Geom_traits::Point_3, Creator>
          g(points[0], points[1], points[2]);
        CGAL::cpp11::copy_n(g, nb_points, out);
      }
    }
  }
  return out;
}

template<class OutputIterator, class TriangleMesh>
OutputIterator
sample_triangle_mesh(const TriangleMesh& tm,
                           double parameter,
                           OutputIterator out,
                           Sampling_method method = RANDOM_UNIFORM)
{
   return  sample_triangle_mesh(
                tm,
                parameter,
                out,
                parameters::all_default(),
                method);
}

template <class Concurrency_tag, class Kernel, class PointRange, class TriangleMesh, class VertexPointMap>
double approximate_Hausdorff_distance(
  const PointRange& original_sample_points,
  const TriangleMesh& tm,
  VertexPointMap vpm
  )
{
  CGAL_assertion_code(  bool is_triangle = is_triangle_mesh(tm) );
  CGAL_assertion_msg (is_triangle,
        "Mesh is not triangulated. Distance computing impossible.");
  #ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "Nb sample points " << sample_points.size() << "\n";
  #endif
  std::vector<typename Kernel::Point_3> sample_points
    (boost::begin(original_sample_points), boost::end(original_sample_points) );

  spatial_sort(sample_points.begin(), sample_points.end());

  typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef AABB_traits<Kernel, Primitive> Traits;
  typedef AABB_tree< Traits > Tree;

  Tree tree( faces(tm).first, faces(tm).second, tm);
  tree.accelerate_distance_queries();
  tree.build();
  typename Kernel::Point_3 hint = get(vpm, *vertices(tm).first);
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<Concurrency_tag,Parallel_tag>::value)
  {
    cpp11::atomic<double> distance(0);
    Distance_computation<Tree, typename Kernel::Point_3> f(tree, hint, sample_points, &distance);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, sample_points.size()), f);
    return distance;
  }
  else
#endif
  {
    double hdist = 0;
    BOOST_FOREACH(const typename Kernel::Point_3& pt, sample_points)
    {
      hint = tree.closest_point(pt, hint);
      typename Kernel::Compute_squared_distance_3 squared_distance;
      typename Kernel::FT dist = squared_distance(hint,pt);
      double d = to_double(CGAL::approximate_sqrt(dist));
      if(d>hdist)
        hdist=d;
    }
    return hdist;
  }
}

template <class Concurrency_tag, class Kernel, class TriangleMesh,
          class NamedParameters,
          class VertexPointMap >
double approximate_Hausdorff_distance(
   const TriangleMesh& m1,
   const TriangleMesh& m2,
   double precision,
   NamedParameters np,
   VertexPointMap vpm,
   Sampling_method method = RANDOM_UNIFORM
)
{
    std::vector<typename Kernel::Point_3> sample_points;
    sample_triangle_mesh(
                m1,
                precision,
                std::back_inserter(sample_points),
                np,
                method );
    return approximate_Hausdorff_distance<Concurrency_tag, Kernel>(sample_points, m2, vpm);
}

// documented functions
/** \ingroup PMP_distance_grp
 * generates points taken on `f`, facet of `tm`, in a way depending on `method` and
 * outputs them to `out`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam OutputIterator a model of `OutputIterator`
 *  holding objects of the same point type as
 *  the value type of the internal vertex point map of `tm`
 * @tparam Sampling_method defines the sampling method. Possible values are:
     - `RANDOM_UNIFORM`: points are generated in a random and uniform way, depending on the area of the face.
     - `GRID`: points are generated on a grid, with a minimum of one point .
     - `MONTE_CARLO`: points are generated randomly . The number of points is proportional to the face area with a minimum of 1.
 *
 * @param tm the triangle mesh that contains `f`
 * @param parameter depends on `method`:
     - `RANDOM_UNIFORM` and `MONTE_CARLO`: the number of points per squared area unit
     - `GRID`: the distance between two consecutive points in the grid
 *
 * @param out output iterator to be filled with sampled points
 * @param np an optional sequence of \ref namedparameters among the ones listed below
 * @param method defines the method of sampling
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map with the points associated to the vertices of `tm`. If this parameter is omitted,
 *    an internal property map for `CGAL::vertex_point_t` should be available for `TriangleMesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPDistanceTraits`\cgalParamEnd
 * \cgalNamedParamsEnd
 */
template<class OutputIterator, class TriangleMesh, class NamedParameters>
OutputIterator
sample_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
            const TriangleMesh& tm,
            double parameter,
            OutputIterator out,
            NamedParameters np,
            Sampling_method method = RANDOM_UNIFORM)
{
    typedef typename GetGeomTraits<TriangleMesh,
            NamedParameters>::type Geom_traits;

    typedef typename GetVertexPointMap<TriangleMesh,
            NamedParameters>::const_type Vpm;

    Vpm pmap = choose_param(get_param(np, vertex_point),
                           get_const_property_map(vertex_point, tm));
    typedef Creator_uniform_3<typename Geom_traits::FT,
                              typename Geom_traits::Point_3> Creator;
  switch(method)
  {
  case RANDOM_UNIFORM:
  {
      std::size_t nb_points = std::ceil( to_double(
                                             face_area(f,tm,parameters::geom_traits(Geom_traits())))*parameter);
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
      typename Geom_traits::Point_3 points[3];
      for(int i=0; i<3; ++i)
      {
          points[i] = get(pmap, target(hd, tm));
          hd = next(hd, tm);
      }
      Random_points_in_triangle_3<typename Geom_traits::Point_3, Creator>
              g(points[0], points[1], points[2]);
      CGAL::cpp11::copy_n(g, nb_points, out);
      return out;
  }
    case GRID:
    {
        //create the triangles and store them
        typename Geom_traits::Point_3 points[3];
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
        for(int i=0; i<3; ++i)
        {
          points[i] = get(pmap, target(hd, tm));
          hd = next(hd, tm);
        }

      internal::triangle_grid_sampling<Geom_traits>(typename Geom_traits::Triangle_3(points[0], points[1], points[2]), parameter, out);
      return out;
    }
    case MONTE_CARLO:
    {
        std::size_t nb_points( (std::max)(
                    std::ceil(to_double(face_area(f,tm,parameters::geom_traits(Geom_traits())))*parameter),
                                         1.) );
        //create the triangles and store them
        typename Geom_traits::Point_3 points[3];
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
        for(int i=0; i<3; ++i)
        {
          points[i] = get(pmap, target(hd, tm));
          hd = next(hd, tm);
        }
        Random_points_in_triangle_3<typename Geom_traits::Point_3, Creator>
          g(points[0], points[1], points[2]);
        CGAL::cpp11::copy_n(g, nb_points, out);
    }
  }
  return out;
}

/**
 * \ingroup PMP_distance_grp
 * computes the approximate Hausdorff distance of `tm1` from `tm2` by
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
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters for `tm1`
 * @tparam NamedParameters2 a sequence of \ref namedparameters for `tm2`
 *
 * @param tm1 the triangle mesh that will be sampled
 * @param tm2 the triangle mesh to compute the distance to
 * @param precision the number of points per squared area unit for the random sampling
 * @param np1 optional sequence of \ref namedparameters for `tm1` among the ones listed below
 * @param np2 optional sequence of \ref namedparameters for `tm2` among the ones listed below
 *
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map with the points associated to the vertices of `tm1` (`tm2`). If this parameter is omitted,
 *    an internal property map for `CGAL::vertex_point_t` should be available for `TriangleMesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPDistanceTraits`\cgalParamEnd
 * \cgalNamedParamsEnd
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximate_Hausdorff_distance( const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        double precision,
                                        const NamedParameters1& np1,
                                        const NamedParameters2& np2)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters1>::type Geom_traits;

  return approximate_Hausdorff_distance<Concurrency_tag, Geom_traits>(
    tm1, tm2,
              precision,
              choose_param(get_param(np1, vertex_point),
                           get_const_property_map(vertex_point, tm1)),

              choose_param(get_param(np2, vertex_point),
                           get_const_property_map(vertex_point, tm2))

  );
}

/**
 * \ingroup PMP_distance_grp
 * computes the approximate symmetric Hausdorff distance between `tm1` and `tm2`.
 * It returns the maximum of `approximate_Hausdorff_distance(tm1,tm2)` and
 * `approximate_Hausdorff_distance(tm2,tm1)`.
 *
 * \copydetails CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance()
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximate_symmetric_Hausdorff_distance(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  double precision,
  const NamedParameters1& np1,
  const NamedParameters2& np2)
{
  return (std::max)(
    approximate_Hausdorff_distance<Concurrency_tag>(tm1,tm2,precision,np1,np2),
    approximate_Hausdorff_distance<Concurrency_tag>(tm2,tm1,precision,np2,np1)
  );
}

/**
 * \ingroup PMP_distance_grp
 * computes the approximate Hausdorff distance between `points` and `tm`
 * @tparam PointRange a range of `Point_3`, model of `Range`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 * @param points the point_set of interest
 * @param tm the triangle mesh to compute the distance to
 * @param np an optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map with the points associated to the vertices of `tm`. If this parameter is omitted,
 *    an internal property map for `CGAL::vertex_point_t` should be available for the
      vertices of `tm` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPDistanceTraits`\cgalParamEnd
 * \cgalNamedParamsEnd
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class PointRange,
          class NamedParameters>
double max_distance_to_triangle_mesh(const PointRange& points,
                                     const TriangleMesh& tm,
                                     const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters>::type Geom_traits;

  return approximate_Hausdorff_distance<Concurrency_tag, Geom_traits>
     (points,tm,choose_param(get_param(np, vertex_point),
                             get_const_property_map(vertex_point, tm)));
}

/*!
 *\ingroup PMP_distance_grp
 *  Computes the approximate Hausdorff distance between `tm` and `points`
 *
 * @tparam PointRange a range of `Point_3`, model of `Range`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 * @param tm the triangle mesh to compute the distance to
 * @param points the point_set of interest.
 * @param precision the precision of the approximated value you want.
 * @param np an optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map with the points associated to the vertices of `tm`. If this parameter is omitted,
 *    an internal property map for `CGAL::vertex_point_t` should be available for the
      vertices of `tm` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPDistanceTraits`. \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template< class TriangleMesh,
          class PointRange,
          class NamedParameters>
double approximate_max_distance_to_point_set(const TriangleMesh& tm,
                                 const PointRange& points,
                                 const double precision,
                                 const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters>::type Geom_traits;

  typedef CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<Geom_traits> > Knn;
  typedef typename Knn::Tree Tree;
  Tree tree(points.begin(), points.end());
  CRefiner<Geom_traits> ref;
  BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(tm))
  {
    typename Geom_traits::Point_3 points[3];
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd(halfedge(f,tm));
    for(int i=0; i<3; ++i)
    {        
        points[i] = get(choose_param(get_param(np, vertex_point),
                                     get_const_property_map(vertex_point, tm)),
                        target(hd, tm));
        hd = next(hd, tm);
    }
    ref.add(points[0], points[1], points[2], tree);
  }
  return to_double(ref.refine(precision, tree));
}

// convenience functions with default parameters

template< class Concurrency_tag,
          class TriangleMesh,
          class PointRange>
double max_distance_to_triangle_mesh(const PointRange& points,
                                     const TriangleMesh& tm)
{
   return max_distance_to_triangle_mesh<Concurrency_tag,
           TriangleMesh,
           PointRange>
           (points, tm, parameters::all_default());
}

template< class TriangleMesh,
          class PointRange>
double approximate_max_distance_to_point_set(const TriangleMesh& tm,
                                 const PointRange& points,
                                 const double precision)
{
  return approximate_max_distance_to_point_set(tm, points, precision, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters>
double approximate_Hausdorff_distance( const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        double precision,
                                        const NamedParameters& np)
{
  return approximate_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision, np, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double approximate_Hausdorff_distance( const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        double precision)
{
  return approximate_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision,
    parameters::all_default(),
    parameters::all_default());
}


template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters>
double approximate_symmetric_Hausdorff_distance(const TriangleMesh& tm1,
                                                 const TriangleMesh& tm2,
                                                 double precision,
                                                 const NamedParameters& np)
{
  return approximate_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision, np, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double approximate_symmetric_Hausdorff_distance(const TriangleMesh& tm1,
                                                 const TriangleMesh& tm2,
                                                 double precision)
{
  return approximate_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, precision,
    parameters::all_default(),
    parameters::all_default());
}

}
} // end of namespace CGAL::Polygon_mesh_processing


#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
