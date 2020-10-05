// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno and Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
#define CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/internal/mesh_to_point_set_hausdorff_distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/utility.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/spatial_sort.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <atomic>
#endif // CGAL_LINKED_WITH_TBB

#include <boost/unordered_set.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <class Kernel, class PointOutputIterator>
PointOutputIterator
triangle_grid_sampling(const typename Kernel::Point_3& p0,
                       const typename Kernel::Point_3& p1,
                       const typename Kernel::Point_3& p2,
                       double distance,
                       PointOutputIterator out)
{
  typename Kernel::Compute_squared_distance_3 squared_distance;
  const double d_p0p1 = to_double(approximate_sqrt( squared_distance(p0, p1) ));
  const double d_p0p2 = to_double(approximate_sqrt( squared_distance(p0, p2) ));

  const double n = (std::max)(std::ceil( d_p0p1 / distance ),
                              std::ceil( d_p0p2 / distance ));

  for(double i=1; i<n; ++i)
  {
    for(double j=1; j<n-i; ++j)
    {
      const double c0=(1-(i+j)/n), c1=i/n, c2=j/n;
      *out++ = typename Kernel::Point_3(p0.x()*c0+p1.x()*c1+p2.x()*c2,
                                        p0.y()*c0+p1.y()*c1+p2.y()*c2,
                                        p0.z()*c0+p1.z()*c1+p2.z()*c2);
    }
  }

  return out;
}

#if defined(CGAL_LINKED_WITH_TBB)
template <class AABB_tree, class PointRange>
struct Distance_computation{
  typedef typename PointRange::const_iterator::value_type Point_3;

  const AABB_tree& tree;
  const PointRange& sample_points;
  Point_3 initial_hint;
  double distance;
  //constructor
  Distance_computation(
          const AABB_tree& tree,
          const Point_3& p,
          const PointRange& sample_points)
    : tree(tree)
    , sample_points(sample_points)
    , initial_hint(p)
    , distance(-1)
  {}
  //split constructor
  Distance_computation(Distance_computation& s, tbb::split )
    : tree(s.tree)
    , sample_points(s.sample_points)
    , initial_hint(s.initial_hint)
    , distance(-1)
  {}

  void
  operator()(const tbb::blocked_range<std::size_t>& range)
  {
    Point_3 hint = initial_hint;
    double hdist = 0;
    for( std::size_t i = range.begin(); i != range.end(); ++i)
    {
      hint = tree.closest_point(*(sample_points.begin() + i), hint);
      typename Kernel_traits<Point_3>::Kernel::Compute_squared_distance_3 squared_distance;
      double d = to_double(CGAL::approximate_sqrt( squared_distance(hint,*(sample_points.begin() + i)) ));
      if(d > hdist)
        hdist=d;
    }
    if(hdist > distance)
      distance = hdist;
  }

  void join( Distance_computation& rhs ) {distance = (std::max)(rhs.distance, distance); }
};
#endif

template <class Concurrency_tag,
          class Kernel,
          class PointRange,
          class AABBTree>
double approximate_Hausdorff_distance_impl(
  const PointRange& sample_points,
  const AABBTree& tree,
  typename Kernel::Point_3 hint)
{
#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(boost::is_convertible<Concurrency_tag,Parallel_tag>::value)
  {
    std::atomic<double> distance;
    distance=0;
    Distance_computation<AABBTree, PointRange> f(tree, hint, sample_points);
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, sample_points.size()), f);
    return f.distance;
  }
  else
#endif
  {
    double hdist = 0;
    for(const typename Kernel::Point_3& pt : sample_points)
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

template<typename PointOutputIterator,
         typename GeomTraits,
         typename NamedParameters,
         typename TriangleIterator,
         typename Randomizer,
         typename Creator,
         typename Derived>
struct Triangle_structure_sampler_base
{
  const NamedParameters np;
  GeomTraits gt;
  PointOutputIterator& out;

  Triangle_structure_sampler_base(PointOutputIterator& out,
                                  const NamedParameters& np)
    : np(np), out(out)
  {}

  void sample_points();
  double get_minimum_edge_length();
  template<typename Tr>
  double get_tr_area(const Tr&);

  template<typename Tr>
  std::array<typename GeomTraits::Point_3, 3> get_tr_points(const Tr& tr);

  void ms_edges_sample(const std::size_t& nb_points_per_edge,
                       const std::size_t& nb_pts_l_u);
  void ru_edges_sample();
  void internal_sample_triangles(double, bool, bool);

  Randomizer get_randomizer();
  std::pair<TriangleIterator, TriangleIterator> get_range();
  std::size_t get_points_size();

  void procede()
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));


    bool use_rs = choose_parameter(get_parameter(np, internal_np::random_uniform_sampling), true);
    bool use_gs = choose_parameter(get_parameter(np, internal_np::grid_sampling), false);
    bool use_ms = choose_parameter(get_parameter(np, internal_np::monte_carlo_sampling), false);

    if(use_gs || use_ms)
      if(is_default_parameter(get_parameter(np, internal_np::random_uniform_sampling)))
        use_rs = false;

    bool smpl_vrtcs = choose_parameter(get_parameter(np, internal_np::do_sample_vertices), true);
    bool smpl_dgs = choose_parameter(get_parameter(np, internal_np::do_sample_edges), true);
    bool smpl_fcs = choose_parameter(get_parameter(np, internal_np::do_sample_faces), true);
    double nb_pts_a_u = choose_parameter(get_parameter(np, internal_np::nb_points_per_area_unit), 0.);
    double nb_pts_l_u = choose_parameter(get_parameter(np, internal_np::nb_points_per_distance_unit), 0.);

    // sample vertices
    if(smpl_vrtcs)
      static_cast<Derived*>(this)->sample_points();

    // grid sampling
    if(use_gs)
    {
      double grid_spacing_ = choose_parameter(get_parameter(np, internal_np::grid_spacing), 0.);

      if(grid_spacing_ == 0.)
      {
        // set grid spacing to the shortest edge length
        grid_spacing_ = static_cast<Derived*>(this)->get_minimum_edge_length();
      }

      static_cast<Derived*>(this)->internal_sample_triangles(grid_spacing_, smpl_fcs, smpl_dgs);
    }

    // monte carlo sampling
    if(use_ms)
    {
      double min_sq_edge_length = (std::numeric_limits<double>::max)();

      std::size_t nb_points_per_face =
          choose_parameter(get_parameter(np, internal_np::number_of_points_per_face), 0);

      std::size_t nb_points_per_edge =
          choose_parameter(get_parameter(np, internal_np::number_of_points_per_edge), 0);

      if((nb_points_per_face == 0 && nb_pts_a_u == 0.) ||
         (nb_points_per_edge == 0 && nb_pts_l_u == 0.))
      {
        min_sq_edge_length = static_cast<Derived*>(this)->get_minimum_edge_length();
      }

      // sample faces
      if(smpl_fcs)
      {
        // set default value
        if(nb_points_per_face == 0 && nb_pts_a_u == 0.)
          nb_pts_a_u = 2. / min_sq_edge_length;

        for(const auto& tr : make_range(static_cast<Derived*>(this)->get_range()))
        {
          std::size_t nb_points = nb_points_per_face;
          if(nb_points == 0)
          {
            nb_points = (std::max)(
                  static_cast<std::size_t>(
                    std::ceil(static_cast<Derived*>(this)->get_tr_area(tr))
                    *nb_pts_a_u), std::size_t(1));
          }

          // extract triangle face points
          std::array<typename GeomTraits::Point_3, 3>points = static_cast<Derived*>(this)->get_tr_points(tr);

          Random_points_in_triangle_3<typename GeomTraits::Point_3, Creator> g(points[0], points[1], points[2]);
          out = std::copy_n(g, nb_points, out);
        }
      }

      // sample edges
      if(smpl_dgs)
        static_cast<Derived*>(this)->ms_edges_sample(nb_points_per_edge, nb_pts_l_u);
    }

    // random uniform sampling
    if(use_rs)
    {
      // sample faces
      if(smpl_fcs)
      {
        std::size_t nb_points
            = choose_parameter(get_parameter(np, internal_np::number_of_points_on_faces), 0);

        typename Derived::Randomizer g = static_cast<Derived*>(this)->get_randomizer();
        if(nb_points == 0)
        {
          if(nb_pts_a_u == 0.)
            nb_points = static_cast<Derived*>(this)->get_points_size();
          else
            nb_points = static_cast<std::size_t>(std::ceil(g.sum_of_weights()*nb_pts_a_u));
        }
        out = std::copy_n(g, nb_points, out);
      }

      // sample edges
      if(smpl_dgs)
        static_cast<Derived*>(this)->ru_edges_sample(nb_pts_l_u,nb_pts_a_u);
    }
  }
};

} // namespace internal

template <class Kernel,
          class FaceRange,
          class TriangleMesh,
          class VertexPointMap,
          class PointOutputIterator>
PointOutputIterator
sample_triangles(const FaceRange& triangles,
                 const TriangleMesh& tm,
                 VertexPointMap vpm,
                 double distance,
                 PointOutputIterator out,
                 bool sample_faces,
                 bool sample_edges,
                 bool add_vertices)
{
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<typename GT::edge_descriptor> sampled_edges;
  boost::unordered_set<typename GT::vertex_descriptor> endpoints;

  for(face_descriptor fd : triangles)
  {
    // sample edges but skip endpoints
    halfedge_descriptor hd = halfedge(fd, tm);
    for(int i=0;i<3; ++i)
    {
      if(sample_edges && sampled_edges.insert(edge(hd, tm)).second )
      {
        Point_ref p0 = get(vpm, source(hd, tm));
        Point_ref p1 = get(vpm, target(hd, tm));
        typename Kernel::Compute_squared_distance_3 squared_distance;
        const double d_p0p1 = to_double(approximate_sqrt(squared_distance(p0, p1)));

        const double nb_pts = std::ceil( d_p0p1 / distance );
        const Vector_3 step_vec =  typename Kernel::Construct_scaled_vector_3()(
          typename Kernel::Construct_vector_3()(p0, p1),
          typename Kernel::FT(1)/typename Kernel::FT(nb_pts));
        for(double i=1; i<nb_pts; ++i)
        {
          *out++=typename Kernel::Construct_translated_point_3()(p0,
            typename Kernel::Construct_scaled_vector_3()(step_vec ,
              typename Kernel::FT(i)));
        }
      }

      //add endpoints once
      if(add_vertices && endpoints.insert(target(hd, tm)).second)
        *out++ = get(vpm, target(hd, tm));

      hd = next(hd, tm);
    }

    // sample triangles
    if(sample_faces)
    {
      Point_ref p0 = get(vpm, source(hd, tm));
      Point_ref p1 = get(vpm, target(hd, tm));
      Point_ref p2 = get(vpm, target(next(hd, tm), tm));
      out = internal::triangle_grid_sampling<Kernel>(p0, p1, p2, distance, out);
    }
  }
  return out;
}

namespace internal {

template<typename Mesh,
         typename PointOutputIterator,
         typename GeomTraits,
         typename Creator,
         typename Vpm,
         typename NamedParameters>
struct Triangle_structure_sampler_for_triangle_mesh
    : Triangle_structure_sampler_base<PointOutputIterator,
                                      GeomTraits,
                                      NamedParameters,
                                      typename boost::graph_traits<Mesh>::face_iterator,
                                      Random_points_in_triangle_mesh_3<Mesh, Vpm, Creator>,
                                      Creator,
                                      Triangle_structure_sampler_for_triangle_mesh<Mesh,
                                                                                   PointOutputIterator,
                                                                                   GeomTraits,
                                                                                   Creator,
                                                                                   Vpm,
                                                                                   NamedParameters> >
{
  typedef Triangle_structure_sampler_for_triangle_mesh<Mesh,
                                                       PointOutputIterator,
                                                       GeomTraits,
                                                       Creator, Vpm,
                                                       NamedParameters>     Self;
  typedef Triangle_structure_sampler_base<PointOutputIterator,
                                          GeomTraits,
                                          NamedParameters,
                                          typename boost::graph_traits<Mesh>::face_iterator,
                                          Random_points_in_triangle_mesh_3<Mesh, Vpm, Creator>,
                                          Creator,
                                          Self>                             Base;

  typedef boost::graph_traits<Mesh>                                         GT;
  typedef typename GT::halfedge_descriptor                                  halfedge_descriptor;
  typedef typename GT::edge_descriptor                                      edge_descriptor;
  typedef typename GT::face_descriptor                                      face_descriptor;

  typedef Random_points_in_triangle_mesh_3<Mesh, Vpm,Creator>               Randomizer;
  typedef typename boost::graph_traits<Mesh>::face_iterator                 TriangleIterator;

  Vpm pmap;
  double min_sq_edge_length;
  const Mesh& tm;

  Triangle_structure_sampler_for_triangle_mesh(const Mesh& m,
                                               PointOutputIterator& out,
                                               const NamedParameters& np)
    : Base(out, np), tm(m)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    pmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                            get_const_property_map(vertex_point, tm));
    min_sq_edge_length = (std::numeric_limits<double>::max)();
  }

  std::pair<TriangleIterator, TriangleIterator> get_range()
  {
    return std::make_pair(faces(tm).begin(), faces(tm).end());
  }

  void sample_points()
  {
    Property_map_to_unary_function<Vpm> unary(pmap);
    this->out = std::copy(boost::make_transform_iterator(boost::begin(vertices(tm)), unary),
                          boost::make_transform_iterator(boost::end(vertices(tm)), unary),
                          this->out);
  }

  double get_minimum_edge_length()
  {
    typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

    if(min_sq_edge_length != (std::numeric_limits<double>::max)())
      return min_sq_edge_length;

    for(edge_descriptor ed : edges(tm))
    {
      const double sq_el = CGAL::to_double(
                             typename GeomTraits::Compute_squared_distance_3()(get(pmap, source(ed, tm)),
                                                                               get(pmap, target(ed, tm))));

      if(sq_el > 0. && sq_el < min_sq_edge_length)
        min_sq_edge_length = sq_el;
    }

    return min_sq_edge_length;
  }

  double get_tr_area(const typename boost::graph_traits<Mesh>::face_descriptor& tr)
  {
    return to_double(face_area(tr,tm,parameters::geom_traits(this->gt)));
  }

  template<typename Tr>//tr = face_descriptor here
  std::array<typename GeomTraits::Point_3, 3> get_tr_points(const Tr& tr)
  {
    std::array<typename GeomTraits::Point_3, 3> points;
    halfedge_descriptor hd(halfedge(tr,tm));
    for(int i=0; i<3; ++i)
    {
      points[i] = get(pmap, target(hd, tm));
      hd = next(hd, tm);
    }
    return points;
  }

  void ms_edges_sample(std::size_t nb_points_per_edge,
                       double nb_pts_l_u)
  {
    typename GeomTraits::Compute_squared_distance_3 squared_distance = this->gt.compute_squared_distance_3_object();

    if(nb_points_per_edge == 0 && nb_pts_l_u == 0.)
      nb_pts_l_u = 1. / CGAL::sqrt(min_sq_edge_length);

    for(edge_descriptor ed : edges(tm))
    {
      std::size_t nb_points = nb_points_per_edge;
      if(nb_points == 0)
      {
        nb_points = (std::max)(
          static_cast<std::size_t>(std::ceil(std::sqrt(to_double(
           squared_distance(get(pmap, source(ed, tm)),
                            get(pmap, target(ed, tm))))) * nb_pts_l_u)),
          std::size_t(1));
      }

      // now do the sampling of the edge
      Random_points_on_segment_3<typename GeomTraits::Point_3, Creator>
        g(get(pmap, source(ed,tm)), get(pmap, target(ed, tm)));
      this->out = std::copy_n(g, nb_points, this->out);
    }
  }
  void ru_edges_sample(double nb_pts_l_u,
                       double nb_pts_a_u)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    std::size_t nb_points = choose_parameter(get_parameter(this->np, internal_np::number_of_points_on_edges), 0);
    Random_points_on_edge_list_graph_3<Mesh, Vpm, Creator> g(tm, pmap);
    if(nb_points == 0)
    {
      if(nb_pts_l_u == 0)
        nb_points = num_vertices(tm);
      else
        nb_points = static_cast<std::size_t>(std::ceil(g.mesh_length() * nb_pts_a_u));
    }
    this->out = std::copy_n(g, nb_points, this->out);
  }

  Randomizer get_randomizer()
  {
    return Randomizer(tm, pmap);
  }

  void internal_sample_triangles(double grid_spacing_, bool smpl_fcs, bool smpl_dgs)
  {
    this->out = sample_triangles<GeomTraits>(faces(tm), tm, pmap, grid_spacing_, this->out, smpl_fcs, smpl_dgs, false);
  }

  std::size_t get_points_size()
  {
    return num_vertices(tm);
  }
};

template<typename PointRange,
         typename TriangleRange,
         typename PointOutputIterator,
         typename GeomTraits,
         typename Creator,
         typename NamedParameters>
struct Triangle_structure_sampler_for_triangle_soup
    : Triangle_structure_sampler_base<PointOutputIterator,
                                      GeomTraits,
                                      NamedParameters,
                                      typename TriangleRange::const_iterator,
                                      Random_points_in_triangle_soup<PointRange,
                                                                     typename TriangleRange::value_type,
                                                                     Creator>,
                                      Creator,
                                      Triangle_structure_sampler_for_triangle_soup<PointRange,
                                                                                   TriangleRange,
                                                                                   PointOutputIterator,
                                                                                   GeomTraits,
                                                                                   Creator,
                                                                                   NamedParameters> >
{
  typedef typename TriangleRange::value_type                                TriangleType;
  typedef Triangle_structure_sampler_for_triangle_soup<PointRange,
                                                       TriangleRange,
                                                       PointOutputIterator,
                                                       GeomTraits,
                                                       Creator,
                                                       NamedParameters>     Self;

  typedef Triangle_structure_sampler_base<PointOutputIterator,
                                          GeomTraits,
                                          NamedParameters,
                                          typename TriangleRange::const_iterator,
                                          Random_points_in_triangle_soup<PointRange, TriangleType, Creator>,
                                          Creator,
                                          Self>                             Base;

  typedef typename GeomTraits::Point_3                                      Point_3;

  typedef Random_points_in_triangle_soup<PointRange, TriangleType, Creator> Randomizer;
  typedef typename TriangleRange::const_iterator                            TriangleIterator;

  double min_sq_edge_length;
  const PointRange& points;
  const TriangleRange& triangles;

  Triangle_structure_sampler_for_triangle_soup(const PointRange& pts,
                                               const TriangleRange& trs,
                                               PointOutputIterator& out,
                                               const NamedParameters& np)
    : Base(out, np), points(pts), triangles(trs)
  {
    min_sq_edge_length = (std::numeric_limits<double>::max)();
  }

  std::pair<TriangleIterator, TriangleIterator> get_range()
  {
    return std::make_pair(triangles.begin(), triangles.end());
  }

  void sample_points()
  {
    this->out = std::copy(points.begin(), points.end(), this->out);
  }

  double get_minimum_edge_length()
  {
    if(min_sq_edge_length != (std::numeric_limits<double>::max)())
      return min_sq_edge_length;

    for(const auto& tr : triangles)
    {
      for(std::size_t i = 0; i< 3; ++i)
      {
        const Point_3& a = points[tr[i]];
        const Point_3& b = points[tr[(i+1)%3]];

        const double sq_el = CGAL::to_double(typename GeomTraits::Compute_squared_distance_3()(a, b));
        if(sq_el > 0. && sq_el < min_sq_edge_length)
          min_sq_edge_length = sq_el;
      }
    }

    return min_sq_edge_length;
  }

  template<typename Tr>
  double get_tr_area(const Tr& tr)
  {
    return to_double(approximate_sqrt(
                       this->gt.compute_squared_area_3_object()(
                         points[tr[0]], points[tr[1]], points[tr[2]])));
  }

  template<typename Tr>
  std::array<Point_3, 3> get_tr_points(const Tr& tr)
  {
    std::array<Point_3, 3> points;
    for(int i=0; i<3; ++i)
    {
      points[i] = this->points[tr[i]];
    }
    return points;
  }

  void ms_edges_sample(std::size_t, double)
  {
    // don't sample edges in soup.
  }

  void ru_edges_sample(double, double)
  {
    // don't sample edges in soup.
  }

  Randomizer get_randomizer()
  {
    return Randomizer(triangles, points);
  }

  void internal_sample_triangles(double distance, bool, bool)
  {
    for(const auto& tr : triangles)
    {
      const Point_3& p0 = points[tr[0]];
      const Point_3& p1 = points[tr[1]];
      const Point_3& p2 = points[tr[2]];

      this->out = internal::triangle_grid_sampling<GeomTraits>(p0, p1, p2, distance, this->out);
    }
  }

  std::size_t get_points_size()
  {
    return points.size();
  }
};

} // namespace internal

/** \ingroup PMP_distance_grp
 *
 * generates points on `tm` and outputs them to `out`; the sampling method
 * is selected using named parameters.
 *
 * @tparam TriangleMesh a model of the concepts `EdgeListGraph` and `FaceListGraph`
 * @tparam PointOutputIterator a model of `OutputIterator`
 *  holding objects of the same point type as
 *  the value type of the point type associated to the mesh `tm`, i.e. the value type of the vertex
 *  point map property map, if provided, or the value type of the internal point property map otherwise
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm the triangle mesh to be sampled
 * @param out output iterator to be filled with sample points
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPDistanceTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_random_uniform_sampling}
 *     \cgalParamDescription{If `true` is passed, points are generated in a random and uniform way
 *                           on the surface of `tm`, and/or on edges of `tm`.}
 *     \cgalParamType{Boolean}
 *     \cgalParamType{`true`}
 *     \cgalParamExtra{For faces, the number of sample points is the value passed to the named
 *                     parameter `number_of_points_on_faces`. If not set,
 *                     the value passed to the named parameter `number_of_points_per_area_unit`
 *                     is multiplied by the area of `tm` to get the number of sample points.
 *                     If none of these parameters is set, the number of points sampled is `num_vertices(tm)`.
 *                     For edges, the number of the number of sample points is the value passed to the named
 *                     parameter `number_of_points_on_edges`. If not set,
 *                     the value passed to the named parameter `number_of_points_per_distance_unit`
 *                     is multiplied by the sum of the length of edges of `tm` to get the number of sample points.
 *                     If none of these parameters is set, the number of points sampled is `num_vertices(tm)`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_grid_sampling}
 *     \cgalParamDescription{If `true` is passed, points are generated on a grid in each triangle,
 *                           with a minimum of one point per triangle.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{The distance between two consecutive points in the grid is that of the length
 *                     of the smallest non-null edge of `tm` or the value passed to the named parameter
 *                     `grid_spacing`. Edges are also split using the same distance, if requested.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_monte_carlo_sampling}
 *     \cgalParamDescription{if `true` is passed, points are generated randomly in each triangle and/or on each edge.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{For faces, the number of points per triangle is the value passed to the named
 *                     parameter `number_of_points_per_face`. If not set, the value passed
 *                     to the named parameter `number_of_points_per_area_unit` is
 *                     used to pick a number of points per face proportional to the triangle
 *                     area with a minimum of one point per face. If none of these parameters
 *                     is set, 2 divided by the square of the length of the smallest non-null
 *                     edge of `tm` is used as if it was passed to
 *                     `number_of_points_per_area_unit`.
 *                     For edges, the number of points per edge is the value passed to the named
 *                     parameter `number_of_points_per_edge`. If not set, the value passed
 *                     to the named parameter `number_of_points_per_distance_unit` is
 *                     used to pick a number of points per edge proportional to the length of
 *                     the edge with a minimum of one point per face. If none of these parameters
 *                     is set, 1 divided by the length of the smallest non-null edge of `tm`
 *                     is used as if it was passed to `number_of_points_per_distance_unit`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{sample_vertices}
 *     \cgalParamDescription{If `true` is passed, the vertices of `tm` are part of the sample.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{sample_edges}
 *     \cgalParamDescription{If `true` is passed, edges of `tm` are sampled.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{sample_faces}
 *     \cgalParamDescription{If `true` is passed, faces of `tm` are sampled.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{grid_spacing}
 *     \cgalParamDescription{a value used as the grid spacing for the grid sampling method}
 *     \cgalParamType{double}
 *     \cgalParamDefault{the length of the shortest, non-degenerate edge of `tm`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_on_edges}
 *     \cgalParamDescription{a value used for the random sampling method as the number of points to pick exclusively on edges}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`num_vertices(tm)` or a value based on `nb_points_per_distance_unit`, if it is defined}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_on_faces}
 *     \cgalParamDescription{a value used for the random sampling method as the number of points to pick on the surface}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`num_vertices(tm)` or a value based on `nb_points_per_area_unit`, if it is defined}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_distance_unit}
 *     \cgalParamDescription{a value used for the random sampling and the Monte Carlo sampling methods to
 *                           respectively determine the total number of points on edges and the number of points per edge}
 *     \cgalParamType{double}
 *     \cgalParamDefault{`1` divided by the length of the shortest, non-degenerate edge of `tm`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_edge}
 *     \cgalParamDescription{a value used by the Monte-Carlo sampling method as the number of points per edge to pick}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`0`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_area_unit}
 *     \cgalParamDescription{a value used for the random sampling and the Monte Carlo sampling methods to
 *                           respectively determine the total number of points inside faces and the number of points per face}
 *     \cgalParamType{double}
 *     \cgalParamDefault{`2` divided by the squared length of the shortest, non-degenerate edge of `tm`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_face}
 *     \cgalParamDescription{a value used by the Monte-Carlo sampling method as the number of points per face to pick}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`0`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @see `CGAL::Polygon_mesh_processing::sample_triangle_soup()`
 */
template<class PointOutputIterator, class TriangleMesh, class NamedParameters>
PointOutputIterator
sample_triangle_mesh(const TriangleMesh& tm,
                     PointOutputIterator out,
                     const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   Vpm;

  internal::Triangle_structure_sampler_for_triangle_mesh<TriangleMesh,
      PointOutputIterator,
      GeomTraits,
      Creator_uniform_3<typename GeomTraits::FT,
      typename GeomTraits::Point_3>,
      Vpm,
      NamedParameters> performer(tm, out, np);
  performer.procede();

  return performer.out;
}

/** \ingroup PMP_distance_grp
 *
 * generates points on a triangle soup and puts them to `out`; the sampling method
 * is selected using named parameters.
 *
 * @tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * @tparam TriangleRange a model of the concept `RandomAccessContainer`
 *                      whose value_type is itself a model of the concept `RandomAccessContainer`
 *                      whose value_type is an unsigned integral value.
 * @tparam PointOutputIterator a model of `OutputIterator` holding objects of the same type as `PointRange`'s value type
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param points the points of the soup
 * @param triangles a `TriangleRange` containing the triangles of the soup to be sampled
 * @param out output iterator to be filled with sample points
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPDistanceTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the point range's point type.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_random_uniform_sampling}
 *     \cgalParamDescription{If `true` is passed, points are generated in a random and uniform way
 *                           over the triangles of the soup.}
 *     \cgalParamType{Boolean}
 *     \cgalParamType{`true`}
 *     \cgalParamExtra{The number of sample points is the value passed to the named
 *                     parameter `number_of_points_on_faces`. If not set,
 *                     the value passed to the named parameter `number_of_points_per_area_unit`
 *                     is multiplied by the area of the soup to get the number of sample points.
 *                     If none of these parameters is set, the number of points sampled is `points.size()`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_grid_sampling}
 *     \cgalParamDescription{If `true` is passed, points are generated on a grid in each triangle,
 *                           with a minimum of one point per triangle.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{The distance between two consecutive points in the grid is that of the length
 *                     of the smallest non-null edge of the soup or the value passed to the named parameter
 *                     `grid_spacing`.}
 *   \cgalParamNEnd
 * *   \cgalParamNBegin{use_monte_carlo_sampling}
 *     \cgalParamDescription{if `true` is passed, points are generated randomly in each triangle.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{The number of points per triangle is the value passed to the named
 *                     parameter `number_of_points_per_face`. If not set, the value passed
 *                     to the named parameter `number_of_points_per_area_unit` is
 *                     used to pick a number of points per face proportional to the triangle
 *                     area with a minimum of one point per face. If none of these parameters
 *                     is set, the number of points per area unit is set to 2 divided
 *                     by the square of the length of the smallest non-null edge of the soup.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{sample_vertices}
 *     \cgalParamDescription{If `true` is passed, the points of `points` are part of the sample.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{sample_faces}
 *     \cgalParamDescription{If `true` is passed, faces of the soup are sampled.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{grid_spacing}
 *     \cgalParamDescription{a value used as the grid spacing for the grid sampling method}
 *     \cgalParamType{double}
 *     \cgalParamDefault{the length of the shortest, non-degenerate edge of the soup}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_on_faces}
 *     \cgalParamDescription{a value used for the random sampling method as the number of points to pick on the surface}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`points.size()` or a value based on `nb_points_per_area_unit`, if it is defined}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_face}
 *     \cgalParamDescription{a value used by the Monte-Carlo sampling method as the number of points per face to pick}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`0`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_points_per_area_unit}
 *     \cgalParamDescription{a value used for the random sampling and the Monte Carlo sampling methods to
 *                           respectively determine the total number of points inside faces and the number of points per face}
 *     \cgalParamType{double}
 *     \cgalParamDefault{`2` divided by the squared length of the shortest, non-degenerate edge of the soup}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \attention Contrary to `sample_triangle_mesh()`, this method does not allow to sample edges.
 *
 * @see `CGAL::Polygon_mesh_processing::sample_triangle_mesh()`
 */
template<class PointOutputIterator,
         class TriangleRange,
         class PointRange,
         class NamedParameters>
PointOutputIterator
sample_triangle_soup(const PointRange& points,
                     const TriangleRange& triangles,
                     PointOutputIterator out,
                     const NamedParameters& np)
{
  typedef typename PointRange::value_type         Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;

  static_assert(std::is_same<Point_3, typename GeomTraits::Point_3>::value, "Wrong point type.");

  internal::Triangle_structure_sampler_for_triangle_soup<PointRange,
      TriangleRange,
      PointOutputIterator,
      GeomTraits,
      Creator_uniform_3<typename GeomTraits::FT,
      typename GeomTraits::Point_3>,
      NamedParameters> performer(points, triangles, out, np);
  performer.procede();

  return performer.out;
}

template<class PointOutputIterator, class TriangleMesh>
PointOutputIterator
sample_triangle_mesh(const TriangleMesh& tm,
                     PointOutputIterator out)
{
  return sample_triangle_mesh(tm, out, parameters::all_default());
}

template<class PointOutputIterator,
         class TriangleRange,
         class PointRange>
PointOutputIterator
sample_triangle_soup(const PointRange& points,
                     const TriangleRange& triangles,
                     PointOutputIterator out)
{
  return sample_triangle_soup(points, triangles, out, parameters::all_default());
}

template <class Concurrency_tag,
          class Kernel,
          class PointRange,
          class TriangleMesh,
          class VertexPointMap>
double approximate_Hausdorff_distance(
  const PointRange& original_sample_points,
  const TriangleMesh& tm,
  VertexPointMap vpm)
{
  CGAL_assertion_code(  bool is_triangle = is_triangle_mesh(tm) );
  CGAL_assertion_msg (is_triangle,
        "Mesh is not triangulated. Distance computing impossible.");
  #ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "Nb sample points " << sample_points.size() << "\n";
  #endif
  typedef typename Kernel::Point_3 Point_3;
  std::vector<Point_3> sample_points
    (boost::begin(original_sample_points), boost::end(original_sample_points) );

  spatial_sort(sample_points.begin(), sample_points.end());

  typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef AABB_tree< AABB_traits<Kernel, Primitive> > Tree;

  Tree tree( faces(tm).first, faces(tm).second, tm);
  tree.build();
  Point_3 hint = get(vpm, *vertices(tm).first);

  return internal::approximate_Hausdorff_distance_impl<Concurrency_tag, Kernel>
    (original_sample_points, tree, hint);
}

template <class Concurrency_tag, class Kernel, class TriangleMesh,
          class NamedParameters,
          class VertexPointMap >
double approximate_Hausdorff_distance(
   const TriangleMesh& tm1,
   const TriangleMesh& tm2,
   const NamedParameters& np,
   VertexPointMap vpm_2)
{
  std::vector<typename Kernel::Point_3> sample_points;
  sample_triangle_mesh(tm1, std::back_inserter(sample_points), np);

  return approximate_Hausdorff_distance<Concurrency_tag, Kernel>(sample_points, tm2, vpm_2);
}

// documented functions

/**
 * \ingroup PMP_distance_grp
 * computes the approximate Hausdorff distance from `tm1` to `tm2` by returning
 * the distance of the farthest point from `tm2` amongst a sampling of `tm1`
 * generated with the function `sample_triangle_mesh()` with
 * `tm1` and `np1` as parameter.
 *
 * A parallel version is provided and requires the executable to be
 * linked against the <a href="https://www.threadingbuildingblocks.org">Intel TBB library</a>.
 * To control the number of threads used, the user may use the `tbb::task_scheduler_init` class.
 * See the <a href="https://www.threadingbuildingblocks.org/documentation">TBB documentation</a>
 * for more details.
 *
 * @tparam Concurrency_tag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam TriangleMesh a model of the concepts `EdgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters" for `tm1`
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters" for `tm2`
 *
 * @param tm1 the triangle mesh that will be sampled
 * @param tm2 the triangle mesh to compute the distance to
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" forwarded to `sample_triangle_mesh()`
 *
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm2`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm2)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * The function `CGAL::parameters::all_default()` can be used to indicate to use the default values
 * for `np1` and specify custom values for `np2`.
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximate_Hausdorff_distance( const TriangleMesh& tm1,
                                       const TriangleMesh& tm2,
                                       const NamedParameters1& np1,
                                       const NamedParameters2& np2)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters1>::type GeomTraits;

  return approximate_Hausdorff_distance<Concurrency_tag, GeomTraits>(
    tm1, tm2, np1, parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_point),
                                get_const_property_map(vertex_point, tm2)));
}

/**
 * \ingroup PMP_distance_grp
 * computes the approximate symmetric Hausdorff distance between `tm1` and `tm2`.
 * It returns the maximum of `approximate_Hausdorff_distance(tm1, tm2, np1, np2)`
 * and `approximate_Hausdorff_distance(tm2, tm1, np2, np1)`.
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double approximate_symmetric_Hausdorff_distance(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const NamedParameters1& np1,
  const NamedParameters2& np2)
{
  return (std::max)(
    approximate_Hausdorff_distance<Concurrency_tag>(tm1,tm2,np1,np2),
    approximate_Hausdorff_distance<Concurrency_tag>(tm2,tm1,np2,np1)
  );
}

/**
 * \ingroup PMP_distance_grp
 * returns the distance to `tm` of the point from `points` that is the furthest from `tm`.
 *
 * @tparam PointRange a range of `Point_3`, model of `Range`. Its iterator type is `RandomAccessIterator`.
 * @tparam TriangleMesh a model of the concepts `EdgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param points the range of points of interest
 * @param tm the triangle mesh to compute the distance to
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPDistanceTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
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
                                 NamedParameters>::type GeomTraits;

  return approximate_Hausdorff_distance<Concurrency_tag, GeomTraits>
     (points,tm,parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tm)));
}

/*!
 *\ingroup PMP_distance_grp
 * returns an approximation of the distance between `points` and the point lying on `tm` that is the farthest from `points`
 *
 * @tparam PointRange a range of `Point_3`, model of `Range`.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm a triangle mesh
 * @param points a range of points
 * @param precision for each triangle of `tm`, the distance of its farthest point from `points` is bounded.
 *                  A triangle is subdivided into sub-triangles so that the difference of its distance bounds
 *                  is smaller than `precision`. `precision` must be strictly positive to avoid infinite loops.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPDistanceTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
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
                                 NamedParameters>::type GeomTraits;
  typedef boost::graph_traits<TriangleMesh> GT;

  typedef Orthogonal_k_neighbor_search<Search_traits_3<GeomTraits> > Knn;
  typedef typename Knn::Tree Tree;
  Tree tree(points.begin(), points.end());
  CRefiner<GeomTraits> ref;
  for(typename GT::face_descriptor f : faces(tm))
  {
    typename GeomTraits::Point_3 points[3];
    typename GT::halfedge_descriptor hd(halfedge(f,tm));
    for(int i=0; i<3; ++i)
    {
      points[i] = get(parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
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
  return approximate_max_distance_to_point_set(tm, points, precision,
                                               parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters>
double approximate_Hausdorff_distance(const TriangleMesh& tm1,
                                      const TriangleMesh& tm2,
                                      const NamedParameters& np)
{
  return approximate_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, np, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double approximate_Hausdorff_distance(const TriangleMesh& tm1,
                                      const TriangleMesh& tm2)
{
  return approximate_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, parameters::all_default(), parameters::all_default());
}


template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters>
double approximate_symmetric_Hausdorff_distance(const TriangleMesh& tm1,
                                                const TriangleMesh& tm2,
                                                const NamedParameters& np)
{
  return approximate_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, np, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double approximate_symmetric_Hausdorff_distance(const TriangleMesh& tm1,
                                                const TriangleMesh& tm2)
{
  return approximate_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1, tm2, parameters::all_default(), parameters::all_default());
}

}
} // end of namespace CGAL::Polygon_mesh_processing


#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
