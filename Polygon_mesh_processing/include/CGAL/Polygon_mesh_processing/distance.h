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
// Author(s)     : Maxime Gimeno, Sebastien Loriot, Martin Skrodzki

#ifndef CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
#define CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/internal/mesh_to_point_set_hausdorff_distance.h>
#include <CGAL/Polygon_mesh_processing/internal/AABB_traversal_traits_with_Hausdorff_distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/utility.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
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
#include <limits>

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
 *                      whose `value_type` is itself a model of the concept `RandomAccessContainer`
 *                      whose `value_type` is an unsigned integral value.
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

////////////////////////////////////////////////////////////////////////

namespace internal {

template <class Concurrency_tag,
          class Kernel,
          class TriangleMesh,
          class VPM1,
          class VPM2>
double bounded_error_Hausdorff_impl(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const typename Kernel::FT& error_bound,
  const bool compare_meshes,
  VPM1 vpm1,
  VPM2 vpm2)
{
  CGAL_assertion_code( const bool is_triangle = is_triangle_mesh(tm1) && is_triangle_mesh(tm2) );
  CGAL_assertion_msg (is_triangle,
        "One of the meshes is not triangulated. Distance computing impossible.");

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM1> TM1_primitive;
  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM2> TM2_primitive;
  typedef AABB_tree< AABB_traits<Kernel, TM1_primitive> > TM1_tree;
  typedef AABB_tree< AABB_traits<Kernel, TM2_primitive> > TM2_tree;
  typedef typename AABB_tree< AABB_traits<Kernel, TM2_primitive> >::AABB_traits Tree_traits;

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  using Candidate = Candidate_triangle<Kernel, face_descriptor>;

  typename Kernel::Compute_squared_distance_3 squared_distance;

  std::vector<face_descriptor> tm1_only, tm2_only;
  std::vector< std::pair<face_descriptor, face_descriptor> > common;

  const auto faces1 = faces(tm1);
  const auto faces2 = faces(tm2);

  CGAL_assertion(faces1.size() > 0);
  CGAL_assertion(faces2.size() > 0);

  TM1_tree tm1_tree;
  TM2_tree tm2_tree;
  if (compare_meshes) {
    match_faces(tm1, tm2, std::back_inserter(common),
      std::back_inserter(tm1_only), std::back_inserter(tm2_only));

    // std::cout << "common t: " << common.size()   << std::endl;
    // std::cout << "tm1 only: " << tm1_only.size() << std::endl;
    // std::cout << "tm2 only: " << tm2_only.size() << std::endl;

    if (tm1_only.size() == 0) return 0.0;
    CGAL_assertion(tm1_only.size() > 0);

    tm1_tree.insert(tm1_only.begin(), tm1_only.end(), tm1, vpm1);
  } else {
    tm1_tree.insert(faces1.begin(), faces1.end(), tm1, vpm1);
  }
  tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);

  // Build an AABB tree on tm1
  tm1_tree.build();
  const bool is_tm1_memory_ok = tm1_tree.accelerate_distance_queries();
  CGAL_assertion(is_tm1_memory_ok);

  // Build an AABB tree on tm2
  tm2_tree.build();
  const bool is_tm2_memory_ok = tm2_tree.accelerate_distance_queries();
  CGAL_assertion(is_tm2_memory_ok);

  // Build traversal traits for tm1_tree
  Hausdorff_primitive_traits_tm1<Tree_traits, Point_3, Kernel, TriangleMesh, VPM1, VPM2> traversal_traits_tm1(
    tm1_tree.traits(), tm2_tree, tm1, tm2, vpm1, vpm2, error_bound );

  // Find candidate triangles in TM1 which might realise the Hausdorff bound
  // TODO Initialize the distances on all the vertices first and store those.
  // TODO Do not traverse TM1, but only TM2, i.e. reduce to Culling on TM2 (Can do this for all triangles in TM1 in parallel)

  tm1_tree.traversal_with_priority( Point_3(0,0,0), traversal_traits_tm1 ); // dummy point given as query as not needed

  // TODO Is there a better/faster data structure than the Heap used here?
  // Can already build a sorted structure while collecting the candidates
  auto& candidate_triangles = traversal_traits_tm1.get_candidate_triangles();
  auto global_bounds = traversal_traits_tm1.get_global_bounds();

  // std::cout << "number of candidate triangles 1: " <<
  //   candidate_triangles.size() << std::endl;

  const FT squared_error_bound = error_bound * error_bound;
  while ( (global_bounds.second - global_bounds.first > error_bound) && !candidate_triangles.empty() ) {

    // Get the first triangle and its Hausdorff bounds from the candidate set
    const Candidate triangle_and_bound = candidate_triangles.top();
    // Remove it from the candidate set as it will be processed now
    candidate_triangles.pop();

    // Only process the triangle if it can contribute to the Hausdorff distance,
    // i.e. if its Upper Bound is higher than the currently known best lower bound
    // and the difference between the bounds to be obtained is larger than the
    // user given error.
    const auto& triangle_bounds = triangle_and_bound.bounds;

    if ( (triangle_bounds.second > global_bounds.first) && (triangle_bounds.second - triangle_bounds.first > error_bound) ) {
      // Get the triangle that is to be subdivided and read its vertices
      const Triangle_3& triangle_for_subdivision = triangle_and_bound.triangle;
      const Point_3 v0 = triangle_for_subdivision.vertex(0);
      const Point_3 v1 = triangle_for_subdivision.vertex(1);
      const Point_3 v2 = triangle_for_subdivision.vertex(2);

      // Check second stopping condition: All three vertices of the triangle
      // are projected onto the same triangle in TM2
      const auto closest_triangle_v0 = tm2_tree.closest_point_and_primitive(v0);
      const auto closest_triangle_v1 = tm2_tree.closest_point_and_primitive(v1);
      const auto closest_triangle_v2 = tm2_tree.closest_point_and_primitive(v2);
      if( (closest_triangle_v0.second == closest_triangle_v1.second) && (closest_triangle_v1.second == closest_triangle_v2.second)) {
        // The upper bound of this triangle is the actual Hausdorff distance of
        // the triangle to the second mesh. Use it as new global lower bound.
        // TODO Update the reference to the realizing triangle here as this is the best current guess.
        global_bounds.first = triangle_bounds.second;
        continue;
      }

      // Check third stopping condition: All edge lengths of the triangle are
      // smaller than the given error bound, cannot get results beyond this
      // bound.
      if (    squared_distance( v0, v1 ) < squared_error_bound
          &&  squared_distance( v0, v2 ) < squared_error_bound
          &&  squared_distance( v1, v2 ) < squared_error_bound ) {
            // The upper bound of this triangle is within error tolerance of
            // the actual upper bound, use it.
            global_bounds.first = triangle_bounds.second;
            continue;
      }

      // Subdivide the triangle into four smaller triangles
      const Point_3 v01 = midpoint( v0, v1 );
      const Point_3 v02 = midpoint( v0, v2 );
      const Point_3 v12 = midpoint( v1, v2 );
      const std::array<Triangle_3,4> sub_triangles = {
        Triangle_3( v0, v01, v02), Triangle_3( v1, v01, v12),
        Triangle_3( v2, v02, v12), Triangle_3( v01, v02, v12)
      };

      // Send each of the four triangles to Culling on B with the bounds of the parent triangle
      for (int i=0; i<4; i++) {
        // Call Culling on B with the single triangle found.
        Hausdorff_primitive_traits_tm2<Tree_traits, Triangle_3, Kernel, TriangleMesh, VPM2> traversal_traits_tm2(
          tm2_tree.traits(), tm2, vpm2,
          triangle_bounds.first,
          triangle_bounds.second,
          infinity_value<FT>(),
          infinity_value<FT>(),
          infinity_value<FT>()
        );
        tm2_tree.traversal_with_priority(sub_triangles[i], traversal_traits_tm2);

        // Update global lower Hausdorff bound according to the obtained local bounds
        const auto local_bounds = traversal_traits_tm2.get_local_bounds();
        if (local_bounds.first > global_bounds.first) {
          global_bounds.first = local_bounds.first;
        }

        // TODO Additionally store the face descriptor of the parent from TM1 in the Candidate_triangle.
        // Add the subtriangle to the candidate list
        candidate_triangles.push(
          Candidate(sub_triangles[i], local_bounds, triangle_and_bound.face) );
      }

      // Update global upper Hausdorff bound after subdivision
      const FT current_max = candidate_triangles.top().bounds.second;
      global_bounds.second = std::max(current_max, global_bounds.first);
    }
  }

  // std::cout << "number of candidate triangles 2: " <<
  //   candidate_triangles.size() << std::endl;

  // Return linear interpolation between found lower and upper bound
  return CGAL::to_double((global_bounds.first + global_bounds.second) / FT(2));

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  // TODO implement parallelized version of the above here.
#endif
}

template <class Point_3,
          class TM2_tree,
          class Kernel>
typename Kernel::FT recursive_hausdorff_subdivision(
  const Point_3& v0,
  const Point_3& v1,
  const Point_3& v2,
  const TM2_tree& tm2_tree,
  const typename Kernel::FT& squared_error_bound)
{
  // If all edge lengths of the triangle are below the error_bound,
  // return maximum of the distances of the three points to TM2 (via TM2_tree).
  const auto max_squared_edge_length =
  std::max(
    std::max(
      squared_distance( v0, v1 ),
      squared_distance( v0, v2 )),
    squared_distance( v1, v2 )
  );
  if ( max_squared_edge_length <  squared_error_bound ) {
    return std::max(
      std::max(
        squared_distance( v0, tm2_tree.closest_point(v0) ),
        squared_distance( v1, tm2_tree.closest_point(v1) ) ),
      squared_distance( v2, tm2_tree.closest_point(v2) )
    );
  }

  // Else subdivide the triangle and proceed recursively
  const Point_3 v01 = midpoint( v0, v1 );
  const Point_3 v02 = midpoint( v0, v2 );
  const Point_3 v12 = midpoint( v1, v2 );

  return std::max (
      std::max(
        recursive_hausdorff_subdivision<Point_3, TM2_tree, Kernel>( v0,v01,v02,tm2_tree,squared_error_bound ),
        recursive_hausdorff_subdivision<Point_3, TM2_tree, Kernel>( v1,v01,v12,tm2_tree,squared_error_bound )
      ),
      std::max(
        recursive_hausdorff_subdivision<Point_3, TM2_tree, Kernel>( v2,v02,v12,tm2_tree,squared_error_bound ),
        recursive_hausdorff_subdivision<Point_3, TM2_tree, Kernel>( v01,v02,v12,tm2_tree,squared_error_bound )
      )
    );
}

template <class Concurrency_tag,
          class Kernel,
          class TriangleMesh,
          class VPM1,
          class VPM2>
double bounded_error_Hausdorff_naive_impl(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const typename Kernel::FT& error_bound,
  VPM1 vpm1,
  VPM2 vpm2)
{
  CGAL_assertion_code( const bool is_triangle = is_triangle_mesh(tm1) && is_triangle_mesh(tm2) );
  CGAL_assertion_msg (is_triangle,
        "One of the meshes is not triangulated. Distance computing impossible.");

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM2> TM2_primitive;
  typedef AABB_tree< AABB_traits<Kernel, TM2_primitive> > TM2_tree;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  // Initially, no lower bound is known
  FT squared_lower_bound = FT(0);
  // Work with squares in the following, only draw sqrt at the very end
  const FT squared_error_bound = error_bound * error_bound;

  // Build an AABB tree on tm2
  TM2_tree tm2_tree( faces(tm2).begin(), faces(tm2).end(), tm2, vpm2 );
  tm2_tree.build();
  tm2_tree.accelerate_distance_queries();

  // Build a map to obtain actual triangles from the face descriptors of tm1.
  const Triangle_from_face_descriptor_map<TriangleMesh, VPM1> face_to_triangle_map( &tm1, vpm1 );

  // Iterate over the triangles of TM1.
  for(const face_descriptor& fd : faces(tm1))
  {
    // Get the vertices of the face and pass them on to a recursive method.
    const Triangle_3 triangle = get(face_to_triangle_map, fd);
    const Point_3 v0 = triangle.vertex(0);
    const Point_3 v1 = triangle.vertex(1);
    const Point_3 v2 = triangle.vertex(2);

    // Recursively process the current triangle to obtain a lower bound on
    // its Hausdorff distance.
    const FT triangle_bound = recursive_hausdorff_subdivision<Point_3, TM2_tree, Kernel>(
      v0, v1, v2, tm2_tree, squared_error_bound );

    // Store the largest lower bound.
    if( triangle_bound > squared_lower_bound ) {
      squared_lower_bound = triangle_bound;
    }
  }

  // Return linear interpolation between found upper and lower bound
  return CGAL::sqrt(CGAL::to_double( squared_lower_bound ));

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  // TODO implement parallelized version of the below here.
#endif
}

} //end of namespace internal

/**
 * \ingroup PMP_distance_grp
 * returns an estimate on the Hausdorff distance between `tm1` and `tm2` that
 * is at most `error_bound` away from the actual Hausdorff distance between
 * the two given meshes.
 * @tparam Concurrency_tag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag`
 *                         and `Parallel_tag`. Currently, parallel computation is
 *                         not implemented, though.
 * @tparam TriangleMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 * @param tm1 a triangle mesh
 * @param tm2 a second triangle mesh
 * @param error_bound Maximum bound by which the Hausdorff distance estimate is
 *                    allowed to deviate from the actual Hausdorff distance.
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` and `tm2` (`np1` and `np2`, respectively)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{match_faces}
 *     \cgalParamDescription{a boolean tag that turns on the preprocessing step that filters out all faces,
 *                           which belong to both meshes and hence do not contribute to the final distance}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *     \cgalParamExtra{Both `np1` and `np2` must have this tag true in order to activate this preprocessing.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double bounded_error_Hausdorff_distance( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound,
                                         const NamedParameters1& np1,
                                         const NamedParameters2& np2)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters1>::type Geom_traits;
  using FT = typename Geom_traits::FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type Vpm1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::const_type Vpm2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  Vpm1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm1));
  Vpm2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm2));

  const bool match_faces1 = parameters::choose_parameter(
    parameters::get_parameter(np1, internal_np::match_faces), true);
  const bool match_faces2 = parameters::choose_parameter(
    parameters::get_parameter(np2, internal_np::match_faces), true);
  const bool match_faces = match_faces1 && match_faces2;

  CGAL_precondition(error_bound >= 0.0);
  const FT error_threshold = static_cast<FT>(error_bound);
  return internal::bounded_error_Hausdorff_impl<Concurrency_tag, Geom_traits>(
    tm1, tm2, error_threshold, match_faces, vpm1, vpm2);
}


template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1>
double bounded_error_Hausdorff_distance( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound,
                                         const NamedParameters1& np1)
{
  return bounded_error_Hausdorff_distance<Concurrency_tag>(tm1, tm2, error_bound, np1, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double bounded_error_Hausdorff_distance( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound)
{
  return bounded_error_Hausdorff_distance<Concurrency_tag>(tm1, tm2, error_bound, parameters::all_default() );
}

/*
 * Implementation of naive Bounded Hausdorff distance computation.
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
double bounded_error_Hausdorff_distance_naive( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound,
                                         const NamedParameters1& np1,
                                         const NamedParameters2& np2)
{
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters1>::type Geom_traits;
  using FT = typename Geom_traits::FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type Vpm1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::const_type Vpm2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  Vpm1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm1));
  Vpm2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm2));

  CGAL_precondition(error_bound >= 0.0);
  const FT error_threshold = static_cast<FT>(error_bound);
  return internal::bounded_error_Hausdorff_naive_impl<Concurrency_tag, Geom_traits>(
    tm1, tm2, error_threshold, vpm1, vpm2);
}

template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1>
double bounded_error_Hausdorff_distance_naive( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound,
                                         const NamedParameters1& np1)
{
  return bounded_error_Hausdorff_distance_naive<Concurrency_tag>(tm1, tm2, error_bound, np1, parameters::all_default());
}

template< class Concurrency_tag,
          class TriangleMesh>
double bounded_error_Hausdorff_distance_naive( const TriangleMesh& tm1,
                                         const TriangleMesh& tm2,
                                         double error_bound)
{
  return bounded_error_Hausdorff_distance_naive<Concurrency_tag>(tm1, tm2, error_bound, parameters::all_default() );
}

} } // end of namespace CGAL::Polygon_mesh_processing


#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
