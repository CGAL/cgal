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
// Author(s)     : Maxime Gimeno, Sebastien Loriot, Martin Skrodzki, Dmitry Anisimov

#ifndef CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
#define CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/internal/mesh_to_point_set_hausdorff_distance.h>
#include <CGAL/Polygon_mesh_processing/internal/AABB_traversal_traits_with_Hausdorff_distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/utility.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Real_timer.h>
#include <CGAL/iterator.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#if defined(CGAL_METIS_ENABLED)
#include <CGAL/boost/graph/partition.h>
#endif // CGAL_METIS_ENABLED

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#endif // CGAL_LINKED_WITH_TBB

#include <any>

#include <unordered_set>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#ifdef CGAL_HAUSDORFF_DEBUG_PP
 #ifndef CGAL_HAUSDORFF_DEBUG
  #define CGAL_HAUSDORFF_DEBUG
 #endif
#endif

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

  const double d_p0p1 = to_double(approximate_sqrt(squared_distance(p0, p1)));
  const double d_p0p2 = to_double(approximate_sqrt(squared_distance(p0, p2)));

  const double n = (std::max)(std::ceil(d_p0p1 / distance),
                              std::ceil(d_p0p2 / distance));

  for(double i=1; i<n; ++i)
  {
    for(double j=1; j<n-i; ++j)
    {
      const double c0=(1-(i+j)/n), c1=i/n, c2=j/n;
      *out++ = typename Kernel::Point_3(p0.x()*c0 + p1.x()*c1 + p2.x()*c2,
                                        p0.y()*c0 + p1.y()*c1 + p2.y()*c2,
                                        p0.z()*c0 + p1.z()*c1 + p2.z()*c2);
    }
  }

  return out;
}

#if defined(CGAL_LINKED_WITH_TBB)
template <class Kernel, class AABB_tree, class PointRange>
struct Distance_computation
{
  typedef typename Kernel::FT FT;
  typedef typename PointRange::const_iterator::value_type Point_3;

  const AABB_tree& tree;
  const PointRange& sample_points;
  Point_3 initial_hint;
  FT sq_distance;

  //constructor
  Distance_computation(const AABB_tree& tree,
                       const Point_3& p,
                       const PointRange& sample_points)
    : tree(tree),
      sample_points(sample_points),
      initial_hint(p),
      sq_distance(-1)
  {}

  //split constructor
  Distance_computation(Distance_computation& s, tbb::split)
    : tree(s.tree),
      sample_points(s.sample_points),
      initial_hint(s.initial_hint),
      sq_distance(-1)
  {}

  void operator()(const tbb::blocked_range<std::size_t>& range)
  {
    Point_3 hint = initial_hint;
    FT sq_hdist = 0;
    typename Kernel_traits<Point_3>::Kernel::Compute_squared_distance_3 squared_distance;

    for(std::size_t i = range.begin(); i != range.end(); ++i)
    {
      hint = tree.closest_point(*(sample_points.begin() + i), hint);
      FT sq_d = squared_distance(hint,*(sample_points.begin() + i));
      if(sq_d > sq_hdist)
        sq_hdist = sq_d;
    }

    if(sq_hdist > sq_distance)
      sq_distance = sq_hdist;
  }

  void join(Distance_computation& rhs) { sq_distance = (std::max)(rhs.sq_distance, sq_distance); }
};
#endif

template <class Concurrency_tag,
          class PointRange,
          class AABBTree,
          class Kernel>
double max_distance_to_mesh_impl(const PointRange& sample_points,
                                 const AABBTree& tree,
                                 typename Kernel::Point_3 hint, // intentional copy
                                 const Kernel& k)
{
  using FT = typename Kernel::FT;

#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!(std::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<Concurrency_tag,Parallel_tag>::value)
  {
    Distance_computation<Kernel, AABBTree, PointRange> f(tree, hint, sample_points);
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, sample_points.size()), f);
    return to_double(approximate_sqrt(f.sq_distance));
  }
  else
#endif
  {
    FT sq_hdist = 0;
    typename Kernel::Compute_squared_distance_3 squared_distance = k.compute_squared_distance_3_object();

    for(const typename Kernel::Point_3& pt : sample_points)
    {
      hint = tree.closest_point(pt, hint);
      FT sq_d = squared_distance(hint, pt);
      if(sq_d > sq_hdist)
        sq_hdist = sq_d;
    }

    return to_double(approximate_sqrt(sq_hdist));
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
  double get_squared_minimum_edge_length();
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
    {
      if(is_default_parameter<NamedParameters, internal_np::random_uniform_sampling_t>::value)
        use_rs = false;
    }

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

      // set grid spacing to the shortest edge length
      if(grid_spacing_ == 0.)
        grid_spacing_ = std::sqrt(static_cast<Derived*>(this)->get_squared_minimum_edge_length());

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
        min_sq_edge_length = static_cast<Derived*>(this)->get_squared_minimum_edge_length();
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
          std::array<typename GeomTraits::Point_3, 3> points = static_cast<Derived*>(this)->get_tr_points(tr);

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

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  std::unordered_set<edge_descriptor> sampled_edges;
  std::unordered_set<vertex_descriptor> endpoints;

  for(face_descriptor fd : triangles)
  {
    // sample edges but skip endpoints
    halfedge_descriptor hd = halfedge(fd, tm);
    for(int i=0;i<3; ++i)
    {
      if(sample_edges && sampled_edges.insert(edge(hd, tm)).second)
      {
        Point_ref p0 = get(vpm, source(hd, tm));
        Point_ref p1 = get(vpm, target(hd, tm));
        typename Kernel::Compute_squared_distance_3 squared_distance;
        const double d_p0p1 = to_double(approximate_sqrt(squared_distance(p0, p1)));

        const double nb_pts = std::ceil(d_p0p1 / distance);
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

  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor           halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor               edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor               face_descriptor;

  typedef typename GeomTraits::FT                                           FT;

  typedef Random_points_in_triangle_mesh_3<Mesh, Vpm,Creator>               Randomizer;
  typedef typename boost::graph_traits<Mesh>::face_iterator                 TriangleIterator;

  Vpm pmap;
  double min_sq_edge_length;
  const Mesh& tm;
  CGAL::Random rnd;

  Triangle_structure_sampler_for_triangle_mesh(const Mesh& m,
                                               PointOutputIterator& out,
                                               const NamedParameters& np)
    : Base(out, np), tm(m)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    CGAL_assertion(!is_empty(tm));

    pmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                            get_const_property_map(vertex_point, tm));

    if(!(is_default_parameter<NamedParameters, internal_np::random_seed_t>::value))
      rnd = CGAL::Random(choose_parameter(get_parameter(np, internal_np::random_seed),0));

    min_sq_edge_length = (std::numeric_limits<double>::max)();
  }

  std::pair<TriangleIterator, TriangleIterator> get_range()
  {
    return std::make_pair(faces(tm).begin(), faces(tm).end());
  }

  void sample_points()
  {
    Property_map_to_unary_function<Vpm> unary(pmap);
    this->out = std::copy(boost::make_transform_iterator(std::begin(vertices(tm)), unary),
                          boost::make_transform_iterator(std::end(vertices(tm)), unary),
                          this->out);
  }

  double get_squared_minimum_edge_length()
  {
    typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

    if(min_sq_edge_length != (std::numeric_limits<double>::max)())
      return min_sq_edge_length;

    FT m_sq_el = min_sq_edge_length;
    for(edge_descriptor ed : edges(tm))
    {
      const FT sq_el = this->gt.compute_squared_distance_3_object()(get(pmap, source(ed, tm)),
                                                                    get(pmap, target(ed, tm)));

      if(sq_el < m_sq_el)
        m_sq_el = sq_el;
    }

    min_sq_edge_length = to_double(m_sq_el);
    return min_sq_edge_length;
  }

  double get_tr_area(const typename boost::graph_traits<Mesh>::face_descriptor& tr)
  {
    return to_double(face_area(tr, tm, parameters::geom_traits(this->gt)));
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
      nb_pts_l_u = 1. / std::sqrt(min_sq_edge_length);

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
    return Randomizer(tm, pmap, rnd);
  }

  void internal_sample_triangles(double grid_spacing_, bool smpl_fcs, bool smpl_dgs)
  {
    this->out = sample_triangles<GeomTraits>(faces(tm), tm, pmap, grid_spacing_,
                                             this->out, smpl_fcs, smpl_dgs, false);
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

  typedef typename GeomTraits::FT                                           FT;
  typedef typename GeomTraits::Point_3                                      Point_3;

  typedef Random_points_in_triangle_soup<PointRange, TriangleType, Creator> Randomizer;
  typedef typename TriangleRange::const_iterator                            TriangleIterator;

  double min_sq_edge_length;
  const PointRange& points;
  const TriangleRange& triangles;
  Random rnd;

  Triangle_structure_sampler_for_triangle_soup(const PointRange& pts,
                                               const TriangleRange& trs,
                                               PointOutputIterator& out,
                                               const NamedParameters& np)
    : Base(out, np), points(pts), triangles(trs)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    min_sq_edge_length = (std::numeric_limits<double>::max)();
    if(!(is_default_parameter<NamedParameters, internal_np::random_seed_t>::value))
      rnd = CGAL::Random(choose_parameter(get_parameter(np, internal_np::random_seed),0));
  }

  std::pair<TriangleIterator, TriangleIterator> get_range()
  {
    return std::make_pair(triangles.begin(), triangles.end());
  }

  void sample_points()
  {
    this->out = std::copy(points.begin(), points.end(), this->out);
  }

  double get_squared_minimum_edge_length()
  {
    if(min_sq_edge_length != (std::numeric_limits<double>::max)())
      return min_sq_edge_length;

    FT m_sq_el = min_sq_edge_length;
    for(const auto& tr : triangles)
    {
      for(std::size_t i = 0; i< 3; ++i)
      {
        const Point_3& a = points[tr[i]];
        const Point_3& b = points[tr[(i+1)%3]];

        const FT sq_el = this->gt.compute_squared_distance_3_object()(a, b);
        if(sq_el < m_sq_el)
          m_sq_el = sq_el;
      }
    }

    min_sq_edge_length = to_double(m_sq_el);
    return min_sq_edge_length;
  }

  template<typename Tr>
  double get_tr_area(const Tr& tr)
  {
    // Kernel_3::Compute_area_3 uses `sqrt()`
    return to_double(approximate_sqrt(
                       this->gt.compute_squared_area_3_object()(
                         points[tr[0]], points[tr[1]], points[tr[2]])));
  }

  template<typename Tr>
  std::array<Point_3, 3> get_tr_points(const Tr& tr)
  {
    std::array<Point_3, 3> points;
    for(int i=0; i<3; ++i)
      points[i] = this->points[tr[i]];

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
    return Randomizer(triangles, points, rnd);
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
 *  the value type of the point type associated to the mesh `tm`, i.e., the value type of the vertex
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
 *   \cgalParamNBegin{random_seed}
 *     \cgalParamDescription{a value to seed the random number generator}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{a value generated with `std::time()`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{use_random_uniform_sampling}
 *     \cgalParamDescription{If `true` is passed, points are generated uniformly at random on faces and/or edges of `tm`.
                             If `do_sample_faces` is `true`, random points will be iteratively generated uniformly at random in the triangle of a face
                             selected with probability proportional to its area. If `do_sample_edges` is `true`, random points will be iteratively generated uniformly at random in the segment of an edge
                             selected with probability proportional to its length.}
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
 *   \cgalParamNBegin{do_sample_vertices}
 *     \cgalParamDescription{If `true` is passed, the vertices of `tm` are part of the sample.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_sample_edges}
 *     \cgalParamDescription{If `true` is passed, edges of `tm` are sampled.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_sample_faces}
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
template<class PointOutputIterator, class TriangleMesh,
         class NamedParameters = parameters::Default_named_parameters>
PointOutputIterator
sample_triangle_mesh(const TriangleMesh& tm,
                     PointOutputIterator out,
                     const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   Vpm;

  CGAL_precondition(!is_empty(tm) && is_triangle_mesh(tm));

  internal::Triangle_structure_sampler_for_triangle_mesh<
      TriangleMesh,
      PointOutputIterator,
      GeomTraits,
      Creator_uniform_3<typename GeomTraits::FT, typename GeomTraits::Point_3>,
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
 *   \cgalParamNBegin{random_seed}
 *     \cgalParamDescription{a value to seed the random number generator}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{a value generated with `std::time()`}
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
 *   \cgalParamNBegin{do_sample_vertices}
 *     \cgalParamDescription{If `true` is passed, the points of `points` are part of the sample.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_sample_faces}
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
         class NamedParameters = parameters::Default_named_parameters>
PointOutputIterator
sample_triangle_soup(const PointRange& points,
                     const TriangleRange& triangles,
                     PointOutputIterator out,
                     const NamedParameters& np = parameters::default_values())
{
  typedef typename PointRange::value_type         Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;

  static_assert(std::is_same<Point_3, typename GeomTraits::Point_3>::value, "Wrong point type.");

  CGAL_precondition(!triangles.empty());

  internal::Triangle_structure_sampler_for_triangle_soup<
      PointRange,
      TriangleRange,
      PointOutputIterator,
      GeomTraits,
      Creator_uniform_3<typename GeomTraits::FT, typename GeomTraits::Point_3>,
      NamedParameters> performer(points, triangles, out, np);
  performer.procede();

  return performer.out;
}

/**
 * \ingroup PMP_distance_grp
 *
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
 *
 * @pre `tm` is a non-empty triangle mesh and `points` is not empty.
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class PointRange,
          class NamedParameters = parameters::Default_named_parameters>
double max_distance_to_triangle_mesh(const PointRange& points,
                                     const TriangleMesh& tm,
                                     const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(!is_empty(tm) && is_triangle_mesh(tm));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GeomTraits;
  typedef typename GeomTraits::Point_3                                            Point_3;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tm));

#ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "Nb sample points " << points.size() << "\n";
#endif

  std::vector<Point_3> points_cpy(std::begin(points), std::end(points));
  spatial_sort(points_cpy.begin(), points_cpy.end());

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef AABB_traits_3<GeomTraits, Primitive> Tree_traits;
  typedef AABB_tree<Tree_traits> Tree;

  Tree_traits tgt/*(gt)*/;
  Tree tree(tgt);
  tree.insert(faces(tm).first, faces(tm).second, tm, vpm);

  const Point_3& hint = get(vpm, *vertices(tm).first);

  return internal::max_distance_to_mesh_impl<Concurrency_tag>(points_cpy, tree, hint, gt);
}

/**
 * \ingroup PMP_distance_grp
 *
 * computes the approximate Hausdorff distance from `tm1` to `tm2` by returning
 * the distance of the farthest point from `tm2` amongst a sampling of `tm1`
 * generated with the function `sample_triangle_mesh()` with
 * `tm1` and `np1` as parameter.
 *
 * A parallel version is provided and requires the executable to be
 * linked against the <a href="https://github.com/oneapi-src/oneTBB">Intel TBB library</a>.
 * To control the number of threads used, the user may use the `tbb::task_scheduler_init` class.
 * See the <a href="https://software.intel.com/content/www/us/en/develop/documentation/onetbb-documentation/top.html">TBB documentation</a>
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
 * @pre `tm1` and `tm2` are non-empty triangle meshes.
 */
template< class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
double approximate_Hausdorff_distance(const TriangleMesh& tm1,
                                      const TriangleMesh& tm2,
                                      const NamedParameters1& np1 = parameters::default_values(),
                                      const NamedParameters2& np2 = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters1>::type GeomTraits;
  typedef typename GeomTraits::Point_3 Point_3;

  CGAL_precondition(!is_empty(tm1) && is_triangle_mesh(tm1));
  CGAL_precondition(!is_empty(tm2) && is_triangle_mesh(tm2));

  std::vector<Point_3> sample_points;
  sample_triangle_mesh(tm1, std::back_inserter(sample_points), np1);

  return max_distance_to_triangle_mesh<Concurrency_tag>(sample_points, tm2, np2);
}

/**
 * \ingroup PMP_distance_grp
 *
 * returns the approximate symmetric Hausdorff distance between `tm1` and `tm2`,
 * that is the maximum of `approximate_Hausdorff_distance(tm1, tm2, np1, np2)`
 * and `approximate_Hausdorff_distance(tm2, tm1, np2, np1)`.
 *
 * See the function `approximate_Hausdorff_distance()` for a complete description of the parameters
 * and requirements.
 */
template <class Concurrency_tag,
          class TriangleMesh,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
double approximate_symmetric_Hausdorff_distance(const TriangleMesh& tm1,
                                                const TriangleMesh& tm2,
                                                const NamedParameters1& np1 = parameters::default_values(),
                                                const NamedParameters2& np2 = parameters::default_values())
{
  return (std::max)(approximate_Hausdorff_distance<Concurrency_tag>(tm1,tm2,np1,np2),
                    approximate_Hausdorff_distance<Concurrency_tag>(tm2,tm1,np2,np1));
}

/*!
 *\ingroup PMP_distance_grp
 *
 * returns an approximation of the distance between `points` and the point lying on `tm` that is the farthest from `points`.
 *
 * @tparam PointRange a range of `Point_3`, model of `Range`
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
 *
 * @pre `tm` is a non-empty triangle mesh and `points` is not empty.
 */
template< class TriangleMesh,
          class PointRange,
          class NamedParameters = parameters::Default_named_parameters>
double approximate_max_distance_to_point_set(const TriangleMesh& tm,
                                             const PointRange& points,
                                             const double precision,
                                             const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(!is_empty(tm) && is_triangle_mesh(tm));
  CGAL_precondition(!points.empty());

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typedef Orthogonal_k_neighbor_search<Search_traits_3<GeomTraits> > Knn;
  typedef typename Knn::Tree Tree;
  Tree tree(points.begin(), points.end());
  CRefiner<GeomTraits> ref;
  for(face_descriptor f : faces(tm))
  {
    typename GeomTraits::Point_3 points[3];
    halfedge_descriptor hd(halfedge(f,tm));
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

////////////////////////////////////////////////////////////////////////

// Use this def in order to get back the parallel version of the one-sided Hausdorff code!
// #define USE_PARALLEL_BEHD

namespace internal {

template <class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2,
          class NamedParameters1,
          class NamedParameters2,
          class TM1Tree,
          class TM2Tree,
          class FaceHandle1,
          class FaceHandle2>
std::pair<typename Kernel::FT, bool>
preprocess_bounded_error_squared_Hausdorff_distance_impl(const TriangleMesh1& tm1,
                                                         const TriangleMesh2& tm2,
                                                         const bool compare_meshes,
                                                         const VPM1 vpm1,
                                                         const VPM2 vpm2,
                                                         const bool is_one_sided_distance,
                                                         const NamedParameters1& np1,
                                                         const NamedParameters2& np2,
                                                         TM1Tree& tm1_tree,
                                                         TM2Tree& tm2_tree,
                                                         std::vector<FaceHandle1>& tm1_only,
                                                         std::vector<FaceHandle2>& tm2_only)
{
  using FT = typename Kernel::FT;

#ifdef CGAL_HAUSDORFF_DEBUG
  using Timer = CGAL::Real_timer;
  Timer timer;
  timer.start();
  std::cout << "* preprocessing begin ...." << std::endl;
  std::cout.precision(17);
#endif

  // Compute the max value that is used as infinity value for the given meshes.
  // In our case, it is twice the length of the diagonal of the bbox of two input meshes.
  const Bbox_3 bbox1 = bbox(tm1);
  const Bbox_3 bbox2 = bbox(tm2);
  const Bbox_3 bb = bbox1 + bbox2;
  const FT sq_dist =   square(bb.xmax() - bb.xmin())
                     + square(bb.ymax() - bb.ymin())
                     + square(bb.zmax() - bb.zmin());

  FT infinity_value = FT(4) * sq_dist;
  CGAL_assertion(infinity_value >= FT(0));

  // Compare meshes and build trees.
  tm1_only.clear();
  tm2_only.clear();
  std::vector<std::pair<FaceHandle1, FaceHandle2> > common;

  const auto faces1 = faces(tm1);
  const auto faces2 = faces(tm2);

  CGAL_precondition(faces1.size() > 0);
  CGAL_precondition(faces2.size() > 0);

  // Compare meshes.
  bool rebuild = false;
  if(compare_meshes) // exact check
  {
    match_faces(tm1, tm2, std::back_inserter(common),
                std::back_inserter(tm1_only), std::back_inserter(tm2_only), np1, np2);

#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout << "-   common: " <<   common.size() << std::endl;
    std::cout << "- tm1 only: " << tm1_only.size() << std::endl;
    std::cout << "- tm2 only: " << tm2_only.size() << std::endl;
#endif

    if(is_one_sided_distance) // one-sided distance
    {
      if(tm1_only.size() > 0) // create TM1 and and full TM2
      {
        tm1_tree.insert(tm1_only.begin(), tm1_only.end(), tm1, vpm1);
        tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);
      }
      else // do not create trees
      {
        CGAL_assertion(tm1_only.size() == 0);
        infinity_value = FT(-1);
      }
    }
    else // symmetric distance
    {
      if(tm1_only.size() == 0 && tm2_only.size() == 0) // do not create trees
      {
        infinity_value = FT(-1);
      }
      else if(common.size() == 0)  // create full TM1 and TM2
      {
        tm1_tree.insert(faces1.begin(), faces1.end(), tm1, vpm1);
        tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);
      }
      else if(tm1_only.size() == 0) // create TM2 and full TM1
      {
        CGAL_assertion(tm2_only.size() > 0);
        CGAL_assertion(tm2_only.size() < faces2.size());
        tm1_tree.insert(faces1.begin(), faces1.end(), tm1, vpm1);
        tm2_tree.insert(tm2_only.begin(), tm2_only.end(), tm2, vpm2);
      }
      else if(tm2_only.size() == 0) // create TM1 and full TM2
      {
        CGAL_assertion(tm1_only.size() > 0);
        CGAL_assertion(tm1_only.size() < faces1.size());
        tm1_tree.insert(tm1_only.begin(), tm1_only.end(), tm1, vpm1);
        tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);
      }
      else // create TM1 and full TM2 and set tag to rebuild them later
      {
        CGAL_assertion(tm1_only.size() > 0);
        CGAL_assertion(tm1_only.size() < faces1.size());
        tm1_tree.insert(tm1_only.begin(), tm1_only.end(), tm1, vpm1);
        tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);
        rebuild = true;
      }
    }
  }
  else // create full TM1 and TM2
  {
    tm1_tree.insert(faces1.begin(), faces1.end(), tm1, vpm1);
    tm2_tree.insert(faces2.begin(), faces2.end(), tm2, vpm2);
  }

#ifdef CGAL_HAUSDORFF_DEBUG
  timer.stop();
  std::cout << "* .... end preprocessing" << std::endl;
  std::cout << "* preprocessing time (sec.): " << timer.time() << std::endl;
#endif

  return std::make_pair(infinity_value, rebuild);
}

template <class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2,
          class TM1Tree,
          class TM2Tree,
          class OutputIterator>
typename Kernel::FT
bounded_error_squared_Hausdorff_distance_impl(const TriangleMesh1& tm1,
                                              const TriangleMesh2& tm2,
                                              const VPM1 vpm1,
                                              const VPM2 vpm2,
                                              const TM1Tree& tm1_tree,
                                              const TM2Tree& tm2_tree,
                                              const typename Kernel::FT error_bound,
                                              const typename Kernel::FT sq_initial_bound,
                                              const typename Kernel::FT sq_distance_bound,
                                              const typename Kernel::FT infinity_value,
                                              OutputIterator& out)
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  auto midpoint = Kernel().construct_midpoint_3_object();

#ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << " -- Bounded Hausdorff --" << std::endl;
  std::cout << "error bound: " << error_bound << std::endl;
  std::cout << "initial bound: " << sq_initial_bound << " (" << approximate_sqrt(sq_initial_bound) << ")" << std::endl;
  std::cout << "distance bound: " << sq_distance_bound << " (" << approximate_sqrt(sq_distance_bound) << ")" << std::endl;
  std::cout << "inf val: " << infinity_value << " (" << approximate_sqrt(infinity_value) << ")" << std::endl;
#endif

  using TM1_hd_traits = Hausdorff_primitive_traits_tm1<Point_3, Kernel, TriangleMesh1, TriangleMesh2, VPM1, VPM2>;
  using TM2_hd_traits = Hausdorff_primitive_traits_tm2<Triangle_3, Kernel, TriangleMesh1, TriangleMesh2, VPM2>;

  using Face_handle_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using Face_handle_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;

  using Candidate = Candidate_triangle<Kernel, Face_handle_1, Face_handle_2>;

  CGAL_precondition(sq_initial_bound >= square(FT(error_bound)));
  CGAL_precondition(sq_distance_bound != FT(0)); // value is -1 if unused
  CGAL_precondition(tm1_tree.size() > 0);
  CGAL_precondition(tm2_tree.size() > 0);

  // First, we apply culling.
#ifdef CGAL_HAUSDORFF_DEBUG
  using Timer = CGAL::Real_timer;
  Timer timer;
  timer.start();
  std::cout << "- applying culling" << std::endl;
  std::cout.precision(17);
#endif

  // Build traversal traits for tm1_tree.
  TM1_hd_traits traversal_traits_tm1(tm2_tree, tm1, tm2, vpm1, vpm2,
                                     infinity_value, sq_initial_bound, sq_distance_bound);

  // Find candidate triangles in TM1, which might realize the Hausdorff bound.
  // We build a sorted structure while collecting the candidates.
  const Point_3 stub(0, 0, 0); // dummy point given as query since it is not needed

  tm1_tree.traversal_with_priority(stub, traversal_traits_tm1);
  auto& candidate_triangles = traversal_traits_tm1.get_candidate_triangles();
  Global_bounds<Kernel, Face_handle_1, Face_handle_2> global_bounds = traversal_traits_tm1.get_global_bounds();

#ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "- bounds post traversal: " << global_bounds.lower << " " << global_bounds.upper << std::endl;
  std::cout << "- number of candidate triangles: " << candidate_triangles.size() << std::endl;
  const FT culling_rate = FT(100) - (FT(candidate_triangles.size()) / FT(tm1_tree.size()) * FT(100));
  std::cout << "- culling rate: " << culling_rate << "%" << std::endl;

  timer.stop();
  std::cout << "* culling (sec.): " << timer.time() << std::endl;
#endif

  CGAL_assertion(global_bounds.lower >= FT(0));
  CGAL_assertion(global_bounds.upper >= global_bounds.lower);
  CGAL_assertion(global_bounds.lpair.first != boost::graph_traits<TriangleMesh1>::null_face());
  CGAL_assertion(global_bounds.lpair.second != boost::graph_traits<TriangleMesh2>::null_face());
  CGAL_assertion(global_bounds.upair.first != boost::graph_traits<TriangleMesh1>::null_face());
  CGAL_assertion(global_bounds.upair.second != boost::graph_traits<TriangleMesh2>::null_face());

  // If we already reached the user-defined max distance bound, we quit.
  if(traversal_traits_tm1.early_exit())
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout << "Quitting early (TM1 traversal): temporary distance " << global_bounds.lower
              << " is already greater than user-defined bound " << sq_distance_bound << std::endl;
#endif

    CGAL_assertion(global_bounds.lower > sq_distance_bound);
    return global_bounds.lower;
  }

  // Second, we apply subdivision.
#ifdef CGAL_HAUSDORFF_DEBUG
  timer.reset();
  std::cout << "- applying subdivision" << std::endl;
  timer.start();
  std::size_t explored_candidates_count = 0;
#endif

  // See Section 5.1 in the paper.
  while(!candidate_triangles.empty())
  {
#ifdef CGAL_HAUSDORFF_DEBUG_PP
    std::cout << "===" << std::endl;
    std::cout << candidate_triangles.size() << " candidates" << std::endl;
    std::cout << "- infinity_value: " << infinity_value << std::endl;
    std::cout << "- error_bound: " << error_bound << std::endl;
    std::cout << "- sq_initial_bound: " << sq_initial_bound << std::endl;
    std::cout << "- sq_distance_bound: " << sq_distance_bound << std::endl;
    std::cout << "- global_bounds.lower: " << global_bounds.lower << std::endl;
    std::cout << "- global_bounds.upper: " << global_bounds.upper << std::endl;
    std::cout << "- diff = " << CGAL::approximate_sqrt(global_bounds.upper) -
                                CGAL::approximate_sqrt(global_bounds.lower) << ", below bound? "
              << ((CGAL::approximate_sqrt(global_bounds.upper) -
                   CGAL::approximate_sqrt(global_bounds.lower)) <= error_bound) << std::endl;
#endif

    CGAL_assertion(global_bounds.lower >= FT(0));
    CGAL_assertion(global_bounds.upper >= global_bounds.lower);

    // @todo could cache those sqrts
    if(CGAL::approximate_sqrt(global_bounds.upper) - CGAL::approximate_sqrt(global_bounds.lower) <= error_bound)
      break;

    // Check if we can early quit.
    if(is_positive(sq_distance_bound)) // empty distance bound is FT(-1)
    {
      const bool early_quit = (sq_distance_bound <= global_bounds.lower);
      if(early_quit)
      {
#ifdef CGAL_HAUSDORFF_DEBUG
        std::cout << "Quitting early with lower bound: " << global_bounds.lower << std::endl;
#endif
        break;
      }
    }

    const Candidate triangle_and_bounds = candidate_triangles.top();
    candidate_triangles.pop();

    // Only process the triangle if it can contribute to the Hausdorff distance,
    // i.e., if its upper bound is higher than the currently known best lower bound
    // and the difference between the bounds to be obtained is larger than the
    // user-given error.
    const auto& triangle_bounds = triangle_and_bounds.bounds;

#ifdef CGAL_HAUSDORFF_DEBUG_PP
    std::cout << "Candidate:" << std::endl;
    std::cout << triangle_and_bounds.triangle.vertex(0) << std::endl;
    std::cout << triangle_and_bounds.triangle.vertex(1) << std::endl;
    std::cout << triangle_and_bounds.triangle.vertex(2) << std::endl;
    std::cout << "triangle_bounds.lower: " << triangle_bounds.lower << std::endl;
    std::cout << "triangle_bounds.upper: " << triangle_bounds.upper << std::endl;
    std::cout << "- diff = " << CGAL::approximate_sqrt(triangle_bounds.upper) -
                                CGAL::approximate_sqrt(triangle_bounds.lower) << ", below bound? "
              << ((CGAL::approximate_sqrt(triangle_bounds.upper) -
                   CGAL::approximate_sqrt(triangle_bounds.lower)) <= error_bound) << std::endl;
#endif

    CGAL_assertion(triangle_bounds.lower >= FT(0));
    CGAL_assertion(triangle_bounds.upper >= triangle_bounds.lower);

    // @todo implement the enclosing-based end criterion (Section 5.1, optional step for TM1 & TM2 closed)

    // Might have been a good candidate when added to the queue, but rendered useless by later insertions
    if(triangle_bounds.upper < global_bounds.lower)
    {
#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "Upper bound is lower than global.lower" << std::endl;
#endif
      continue;
    }

    if((CGAL::approximate_sqrt(triangle_bounds.upper) - CGAL::approximate_sqrt(triangle_bounds.lower)) <= error_bound)
    {
#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "Candidate triangle bounds are tight enough: " << triangle_bounds.lower << " " << triangle_bounds.upper << std::endl;
#endif
      continue;
    }

#ifdef CGAL_HAUSDORFF_DEBUG
    ++explored_candidates_count;
#endif

    // Triangle to be subdivided
    const Triangle_3& triangle_for_subdivision = triangle_and_bounds.triangle;
    const Point_3& v0 = triangle_for_subdivision.vertex(0);
    const Point_3& v1 = triangle_for_subdivision.vertex(1);
    const Point_3& v2 = triangle_for_subdivision.vertex(2);

    // Stopping condition: All three vertices of the triangle are projected onto the same triangle in TM2.
    const auto closest_triangle_v0 = tm2_tree.closest_point_and_primitive(v0);
    const auto closest_triangle_v1 = tm2_tree.closest_point_and_primitive(v1);
    const auto closest_triangle_v2 = tm2_tree.closest_point_and_primitive(v2);
    CGAL_assertion(closest_triangle_v0.second != boost::graph_traits<TriangleMesh2>::null_face());
    CGAL_assertion(closest_triangle_v1.second != boost::graph_traits<TriangleMesh2>::null_face());
    CGAL_assertion(closest_triangle_v2.second != boost::graph_traits<TriangleMesh2>::null_face());

    if((closest_triangle_v0.second == closest_triangle_v1.second) &&
       (closest_triangle_v1.second == closest_triangle_v2.second))
    {
#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "Projects onto the same TM2 face" << std::endl;
#endif

      // The upper bound of this triangle is the actual Hausdorff distance of
      // the triangle to the second mesh. Use it as new global lower bound.
      // Here, we update the reference to the realizing triangle as this is the best current guess.
      global_bounds.lower = triangle_bounds.upper;
      global_bounds.lpair.second = triangle_bounds.tm2_uface;

      continue;
    }

    // Subdivide the triangle into four smaller triangles.
    const Point_3 v01 = midpoint(v0, v1);
    const Point_3 v02 = midpoint(v0, v2);
    const Point_3 v12 = midpoint(v1, v2);
    const std::array<Triangle_3, 4> sub_triangles = { Triangle_3(v0, v01, v02), Triangle_3(v1 , v01, v12),
                                                      Triangle_3(v2, v02, v12), Triangle_3(v01, v02, v12) };

    // Send each of the four triangles to culling on B
    for(std::size_t i=0; i<4; ++i)
    {
      // Call culling on B with the single triangle found.
#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "\nSubface #" << i << "\n"
                << "Geometry: " << sub_triangles[i] << std::endl;
#endif

      // Checking as in during TM1 culling is expensive

      // @todo? For each sub-triangle `ts1` that has a vertex of `v` of the triangle `t1` being subdivided,
      // we have a lower bound on `h(ts1, TM2)` because:
      //  h_t1_lower = max_{vi in t1} min_{t2 in TM2} d(vi, t2)
      // and
      //  h_ts1_lower = max_{vi in ts1} min_{t2 in TM2} d(vi, t2) > min_{t2 in TM2} d(v, t2)
      // But:
      // - we don't keep that in memory (not very hard to change, simply put `m_hi_lower`
      //   from the TM2 traversal traits into the candidate
      // - what's the point? TM2 culling is performed on the local upper bound, so is there
      //   a benefit from providing this value?
      //
      // (We also have that error_bound is a lower bound.)
      const Bbox_3 sub_t1_bbox = sub_triangles[i].bbox();

      // The lower bound is:
      //   h_lower(t1, TM2) := max_{v in t1} min_{t2 in TM2} d(v, t2)

      // The upper bound is:
      //   h_upper(t1, TM2) := min_{t2 in TM2} max_{v in t1} d(v, t2)
      // The value max_{p in t1} d(p, t2) is realized at a vertex of t1.
      // Thus, when splitting t1 into four subtriangles, the distance at the three new vertices
      // is smaller than max_{v in t1} d(v, t2)
      // Thus, subdivision can only decrease the min, and the upper bound.
      Local_bounds<Kernel, Face_handle_1, Face_handle_2> bounds(triangle_bounds.upper);

      // Ensure 'uface' is initialized in case the upper bound is not changed by the subdivision
      bounds.tm2_uface = triangle_bounds.tm2_uface;

      TM2_hd_traits traversal_traits_tm2(sub_t1_bbox, tm2, vpm2, bounds, global_bounds, infinity_value);
      tm2_tree.traversal_with_priority(sub_triangles[i], traversal_traits_tm2);

      // Update global lower Hausdorff bound according to the obtained local bounds.
      const auto& sub_triangle_bounds = traversal_traits_tm2.get_local_bounds();

#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "Subdivided triangle bounds: " << sub_triangle_bounds.lower << " " << sub_triangle_bounds.upper << std::endl;
#endif

      CGAL_assertion(sub_triangle_bounds.lower >= FT(0));
      CGAL_assertion(sub_triangle_bounds.upper >= sub_triangle_bounds.lower);
      CGAL_assertion(sub_triangle_bounds.tm2_lface != boost::graph_traits<TriangleMesh2>::null_face());
      CGAL_assertion(sub_triangle_bounds.tm2_uface != boost::graph_traits<TriangleMesh2>::null_face());

      // The global lower bound is the max of the per-face lower bounds
      if(sub_triangle_bounds.lower > global_bounds.lower)
      {
        global_bounds.lower = sub_triangle_bounds.lower;
        global_bounds.lpair.first = triangle_and_bounds.tm1_face;
        global_bounds.lpair.second = sub_triangle_bounds.tm2_lface;
      }

      // The global upper bound is:
      //   max_{query in TM1} min_{primitive in TM2} max_{v in query} (d(v, primitive))
      // which can go down, so it is only recomputed once splitting is finished,
      // using the top value of the PQ

      candidate_triangles.emplace(sub_triangles[i], sub_triangle_bounds, triangle_and_bounds.tm1_face);
    }

    // Update global upper Hausdorff bound after subdivision.
    const Candidate& top_candidate = candidate_triangles.top();
    const FT current_upmost = top_candidate.bounds.upper;
#ifdef CGAL_HAUSDORFF_DEBUG_PP
    std::cout << "global_bounds.lower = " << global_bounds.lower << std::endl;
    std::cout << "global_bounds.upper = " << global_bounds.upper << std::endl;
    std::cout << "current upper bound = " << current_upmost << std::endl;
#endif

    CGAL_assertion(is_positive(current_upmost));

    if(current_upmost < global_bounds.lower)
    {
#ifdef CGAL_HAUSDORFF_DEBUG_PP
      std::cout << "Top of the queue is lower than the lowest!" << std::endl;
#endif

      global_bounds.upper = global_bounds.lower; // not really needed since lower is returned but doesn't hurt
      global_bounds.upair.first = global_bounds.lpair.first;
      global_bounds.upair.second = global_bounds.lpair.second;

      break;
    }

    CGAL_assertion(current_upmost >= global_bounds.lower);

    global_bounds.upper = current_upmost;
    global_bounds.upair.first = top_candidate.tm1_face;
    global_bounds.upair.second = top_candidate.bounds.tm2_uface;

#ifdef CGAL_HAUSDORFF_DEBUG_PP
    std::cout << "Global bounds post subdi: " << global_bounds.lower << " " << global_bounds.upper << std::endl;
#endif

    CGAL_assertion(global_bounds.lower >= FT(0));
    CGAL_assertion(global_bounds.upper >= global_bounds.lower);
  }

#ifdef CGAL_HAUSDORFF_DEBUG
  timer.stop();
  std::cout << "* subdivision (sec.): " << timer.time() << std::endl;
  std::cout << "Explored " << explored_candidates_count << " candidates" << std::endl;
  std::cout << "Final global bounds: " << global_bounds.lower << " " << global_bounds.upper << std::endl;
  std::cout << "Final global bounds (sqrt): " << CGAL::approximate_sqrt(global_bounds.lower) << " "
                                              << CGAL::approximate_sqrt(global_bounds.upper) << std::endl;
  std::cout << "Difference: " << CGAL::approximate_sqrt(global_bounds.upper) -
                                 CGAL::approximate_sqrt(global_bounds.lower) << std::endl;
#endif

  CGAL_assertion(global_bounds.lower >= FT(0));
  CGAL_assertion(global_bounds.upper >= global_bounds.lower);
  CGAL_assertion(CGAL::approximate_sqrt(global_bounds.upper) - CGAL::approximate_sqrt(global_bounds.lower) <= error_bound);

  // Get realizing triangles.
  CGAL_assertion(global_bounds.lpair.first != boost::graph_traits<TriangleMesh1>::null_face());
  CGAL_assertion(global_bounds.lpair.second != boost::graph_traits<TriangleMesh2>::null_face());
  CGAL_assertion(global_bounds.upair.first != boost::graph_traits<TriangleMesh1>::null_face());
  CGAL_assertion(global_bounds.upair.second != boost::graph_traits<TriangleMesh2>::null_face());

  // Output face pairs, which realize the Hausdorff distance.
  *out++ = global_bounds.lpair;
  *out++ = global_bounds.upair;

  // Return the lower bound because if the correct value is in [0; lower_bound[, the result
  // must still be within the error bound (we have set lower_bound to error_bound initially)
  return global_bounds.lower;
}

#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)

template<class TriangleMesh, class VPM, class TMTree>
struct Triangle_mesh_wrapper
{
  const TriangleMesh& tm; const VPM& vpm;
  const bool is_tm2; TMTree& tm_tree;
  Triangle_mesh_wrapper(const TriangleMesh& tm, const VPM& vpm,
                        const bool is_tm2, TMTree& tm_tree)
    : tm(tm), vpm(vpm), is_tm2(is_tm2), tm_tree(tm_tree)
  { }

  void build_tree()
  {
    tm_tree.insert(faces(tm).begin(), faces(tm).end(), tm, vpm);
    tm_tree.build();
    if(is_tm2)
      tm_tree.accelerate_distance_queries();
    else
      tm_tree.do_not_accelerate_distance_queries();
  }
};

template<class TM1Wrapper, class TM2Wrapper>
struct Bounded_error_preprocessing
{
#ifdef CGAL_HAUSDORFF_DEBUG
  using Timer = CGAL::Real_timer;
#endif
  std::vector<std::any>& tm_wrappers;

  // Constructor.
  Bounded_error_preprocessing(std::vector<std::any>& tm_wrappers)
    : tm_wrappers(tm_wrappers)
  { }

  // Split constructor.
  Bounded_error_preprocessing(Bounded_error_preprocessing& s, tbb::split)
    : tm_wrappers(s.tm_wrappers)
  { }

  bool is_tm1_wrapper(const std::any& operand) const { return operand.type() == typeid(TM1Wrapper); }
  bool is_tm2_wrapper(const std::any& operand) const { return operand.type() == typeid(TM2Wrapper); }

  // TODO: make AABB tree build parallel!
  void operator()(const tbb::blocked_range<std::size_t>& range)
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    Timer timer;
    timer.reset();
    timer.start();
    std::cout.precision(17);
#endif

    for(std::size_t i = range.begin(); i != range.end(); ++i)
    {
      CGAL_assertion(i < tm_wrappers.size());
      auto& tm_wrapper = tm_wrappers[i];
      if(is_tm1_wrapper(tm_wrapper))
      {
        TM1Wrapper& object = std::any_cast<TM1Wrapper&>(tm_wrapper);
        object.build_tree();
      }
      else if(is_tm2_wrapper(tm_wrapper))
      {
        TM2Wrapper& object = std::any_cast<TM2Wrapper&>(tm_wrapper);
        object.build_tree();
      }
      else
      {
        CGAL_assertion_msg(false, "Error: wrong boost any type!");
      }
    }

#ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    std::cout << "* time operator() preprocessing (sec.): " << timer.time() << std::endl;
#endif
  }

  void join(Bounded_error_preprocessing&) { }
};

template <class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2,
          class TM1Tree,
          class TM2Tree,
          class Kernel>
struct Bounded_error_squared_distance_computation
{
  using FT = typename Kernel::FT;
#ifdef CGAL_HAUSDORFF_DEBUG
  using Timer = CGAL::Real_timer;
#endif

  const std::vector<TriangleMesh1>& tm1_parts;
  const TriangleMesh2& tm2;
  const double error_bound;
  const VPM1 vpm1; const VPM2 vpm2;
  const FT infinity_value;
  const FT sq_initial_bound;
  const std::vector<TM1Tree>& tm1_trees;
  const TM2Tree& tm2_tree;
  FT sq_hdist;

  // Constructor.
  Bounded_error_squared_distance_computation(const std::vector<TriangleMesh1>& tm1_parts,
                                             const TriangleMesh2& tm2,
                                             const double error_bound,
                                             const VPM1 vpm1, const VPM2 vpm2,
                                             const FT infinity_value,
                                             const FT sq_initial_bound,
                                             const std::vector<TM1Tree>& tm1_trees,
                                             const TM2Tree& tm2_tree)
    : tm1_parts(tm1_parts), tm2(tm2),
      error_bound(error_bound),
      vpm1(vpm1), vpm2(vpm2),
      infinity_value(infinity_value), sq_initial_bound(sq_initial_bound),
      tm1_trees(tm1_trees), tm2_tree(tm2_tree),
      sq_hdist(-1)
  {
    CGAL_assertion(tm1_parts.size() == tm1_trees.size());
  }

  // Split constructor.
  Bounded_error_squared_distance_computation(Bounded_error_squared_distance_computation& s, tbb::split)
    : tm1_parts(s.tm1_parts), tm2(s.tm2),
      error_bound(s.error_bound),
      vpm1(s.vpm1), vpm2(s.vpm2),
      infinity_value(s.infinity_value), sq_initial_bound(s.sq_initial_bound),
      tm1_trees(s.tm1_trees), tm2_tree(s.tm2_tree),
      sq_hdist(-1)
  {
    CGAL_assertion(tm1_parts.size() == tm1_trees.size());
  }

  void operator()(const tbb::blocked_range<std::size_t>& range)
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    Timer timer;
    timer.reset();
    timer.start();
    std::cout.precision(17);
#endif

    FT sq_dist = FT(-1);
    auto stub = CGAL::Emptyset_iterator();

    for(std::size_t i = range.begin(); i != range.end(); ++i)
    {
      CGAL_assertion(i < tm1_parts.size());
      CGAL_assertion(i < tm1_trees.size());
      const auto& tm1 = tm1_parts[i];
      const auto& tm1_tree = tm1_trees[i];

      // TODO: add distance_bound (now it is FT(-1)) in case we use parallel
      // for checking if two meshes are close.
      const FT sqd = bounded_error_squared_Hausdorff_distance_impl<Kernel>(
                       tm1, tm2, vpm1, vpm2, tm1_tree, tm2_tree,
                       error_bound, sq_initial_bound, FT(-1) /*sq_distance_bound*/, infinity_value,
                       stub);
      if(sqd > sq_dist)
        sq_dist = sqd;
    }

    if(sq_dist > sq_hdist)
      sq_hdist = sq_dist;

#ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    std::cout << "* time operator() computation (sec.): " << timer.time() << std::endl;
#endif
  }

  void join(Bounded_error_squared_distance_computation& rhs)
  {
    sq_hdist = (CGAL::max)(rhs.sq_hdist, sq_hdist);
  }
};

#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)

template <class Concurrency_tag,
          class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2,
          class NamedParameters1,
          class NamedParameters2,
          class OutputIterator>
typename Kernel::FT
bounded_error_squared_one_sided_Hausdorff_distance_impl(const TriangleMesh1& tm1,
                                                        const TriangleMesh2& tm2,
                                                        const typename Kernel::FT error_bound,
                                                        const typename Kernel::FT sq_distance_bound,
                                                        const bool compare_meshes,
                                                        const VPM1 vpm1,
                                                        const VPM2 vpm2,
                                                        const NamedParameters1& np1,
                                                        const NamedParameters2& np2,
                                                        OutputIterator& out)
{
#if !defined(CGAL_LINKED_WITH_TBB) || !defined(CGAL_METIS_ENABLED)
  static_assert(!std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value,
                                      "Parallel_tag is enabled but at least TBB or METIS is unavailable.");
#endif

  using FT = typename Kernel::FT;

  using TM1 = TriangleMesh1;
  using TM2 = TriangleMesh2;

  using TM1_primitive = AABB_face_graph_triangle_primitive<TM1, VPM1>;
  using TM2_primitive = AABB_face_graph_triangle_primitive<TM2, VPM2>;

  using TM1_traits = AABB_traits_3<Kernel, TM1_primitive>;
  using TM2_traits = AABB_traits_3<Kernel, TM2_primitive>;

  using TM1_tree = AABB_tree<TM1_traits>;
  using TM2_tree = AABB_tree<TM2_traits>;

  using Face_handle_1 = typename boost::graph_traits<TM1>::face_descriptor;
  using Face_handle_2 = typename boost::graph_traits<TM2>::face_descriptor;

  // This is parallel version: we split the tm1 into parts, build trees for all parts, and
  // run in parallel all BHD computations. The final distance is obtained by taking the max
  // between BHDs computed for these parts with respect to tm2.
  // This is off by default because the parallel version does not show much of runtime improvement.
  // The slowest part is building AABB trees and this is what should be accelerated in the future.
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)
  using TMF           = CGAL::Face_filtered_graph<TM1>;
  using TMF_primitive = AABB_face_graph_triangle_primitive<TMF, VPM1>;
  using TMF_traits    = AABB_traits_3<Kernel, TMF_primitive>;
  using TMF_tree      = AABB_tree<TMF_traits>;
  using TM1_wrapper   = Triangle_mesh_wrapper<TMF, VPM1, TMF_tree>;
  using TM2_wrapper   = Triangle_mesh_wrapper<TM2, VPM2, TM2_tree>;

  std::vector<TMF> tm1_parts;
  std::vector<TMF_tree> tm1_trees;
  std::vector<std::any> tm_wrappers;
#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)

#ifdef CGAL_HAUSDORFF_DEBUG
  using Timer = CGAL::Real_timer;
  Timer timer;
  std::cout.precision(17);
#endif

  TM1_tree tm1_tree;
  TM2_tree tm2_tree;

  FT infinity_value = FT(-1);

#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)
  // TODO: add to NP!
  const int nb_cores = 4;
  const std::size_t min_nb_faces_to_split = 100; // TODO: increase this number?
 #ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "* num cores: " << nb_cores << std::endl;
 #endif

  if(std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value &&
     nb_cores > 1 &&
     faces(tm1).size() >= min_nb_faces_to_split)
  {
    // (0) -- Compute infinity value.
 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.reset();
    timer.start();
 #endif

    const Bbox_3 bbox1 = bbox(tm1);
    const Bbox_3 bbox2 = bbox(tm2);
    const Bbox_3 bb = bbox1 + bbox2;
    const FT sq_dist =   square(bb.xmax() - bb.xmin())
                       + square(bb.ymax() - bb.ymin())
                       + square(bb.zmax() - bb.zmin());
    infinity_value = FT(4) * sq_dist;
    CGAL_assertion(infinity_value >= FT(0));

 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    const double time0 = timer.time();
    std::cout << "- computing infinity (sec.): " << time0 << std::endl;
 #endif

    // (1) -- Create partition of tm1.
 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.reset();
    timer.start();
 #endif

    using Face_property_tag = CGAL::dynamic_face_property_t<int>;
    auto face_pid_map = get(Face_property_tag(), tm1);
    CGAL::METIS::partition_graph(tm1, nb_cores, CGAL::parameters::face_partition_id_map(face_pid_map));

 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    const double time1 = timer.time();
    std::cout << "- computing partition time (sec.): " << time1 << std::endl;
 #endif

    // (2) -- Create a filtered face graph for each part.
 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.reset();
    timer.start();
 #endif

    tm1_parts.reserve(nb_cores);
    for(int i = 0; i < nb_cores; ++i)
    {
      tm1_parts.emplace_back(tm1, i, face_pid_map);
      // TODO: why is it triggered sometimes?
      // CGAL_assertion(tm1_parts.back().is_selection_valid());
 #ifdef CGAL_HAUSDORFF_DEBUG
      std::cout << "- part " << i << " size: " << tm1_parts.back().number_of_faces() << std::endl;
 #endif
    }

    CGAL_assertion(tm1_parts.size() == nb_cores);
 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    const double time2 = timer.time();
    std::cout << "- creating graphs time (sec.): " << time2 << std::endl;
 #endif

    // (3) -- Preprocess all input data.
 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.reset();
    timer.start();
 #endif

    tm1_trees.resize(tm1_parts.size());
    tm_wrappers.reserve(tm1_parts.size() + 1);
    for(std::size_t i = 0; i < tm1_parts.size(); ++i)
      tm_wrappers.push_back(TM1_wrapper(tm1_parts[i], vpm1, false, tm1_trees[i]));

    tm_wrappers.push_back(TM2_wrapper(tm2, vpm2, true, tm2_tree));
    CGAL_assertion(tm_wrappers.size() == tm1_parts.size() + 1);

    Bounded_error_preprocessing<TM1_wrapper, TM2_wrapper> bep(tm_wrappers);
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, tm_wrappers.size()), bep);

 #ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    const double time3 = timer.time();
    std::cout << "- creating trees time (sec.) " << time3 << std::endl;
 #endif

 #ifdef CGAL_HAUSDORFF_DEBUG
    // Final timing
    std::cout << "* preprocessing parallel time (sec.) " << time0 + time1 + time2 + time3 << std::endl;
 #endif

  } else // sequential version
#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    timer.reset();
    timer.start();
    std::cout << "* preprocessing sequential version " << std::endl;
#endif

    bool rebuild = false;
    std::vector<Face_handle_1> tm1_only;
    std::vector<Face_handle_2> tm2_only;
    std::tie(infinity_value, rebuild) =
        preprocess_bounded_error_squared_Hausdorff_distance_impl<Kernel>(
          tm1, tm2, compare_meshes, vpm1, vpm2, true /*is_one_sided_distance*/, np1, np2,
          tm1_tree, tm2_tree, tm1_only, tm2_only);

    CGAL_assertion(!rebuild);

    if(infinity_value >= FT(0))
    {
      tm1_tree.build();
      tm2_tree.build();
      tm1_tree.do_not_accelerate_distance_queries();
      tm2_tree.accelerate_distance_queries();
    }

#ifdef CGAL_HAUSDORFF_DEBUG
    timer.stop();
    std::cout << "* preprocessing sequential time (sec.) " << timer.time() << std::endl;
#endif
  }

#ifdef CGAL_HAUSDORFF_DEBUG
  std::cout << "* infinity_value: " << infinity_value << std::endl;
#endif

  if(is_negative(infinity_value))
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout << "* culling rate: 100%" << std::endl;
#endif
    const auto face1 = *(faces(tm1).begin());
    const auto face2 = *(faces(tm2).begin());
    *out++ = std::make_pair(face1, face2);
    *out++ = std::make_pair(face1, face2);
    return 0.; // TM1 is part of TM2 so the distance is zero
  }

  CGAL_assertion(infinity_value > FT(0));
  CGAL_assertion(error_bound >= 0.);

  const FT sq_initial_bound = square(FT(error_bound));
  FT sq_hdist = FT(-1);

#ifdef CGAL_HAUSDORFF_DEBUG
  timer.reset();
  timer.start();
#endif

#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED) && defined(USE_PARALLEL_BEHD)
  if(std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value &&
     nb_cores > 1 &&
     faces(tm1).size() >= min_nb_faces_to_split)
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout << "* executing parallel version " << std::endl;
#endif

    using Comp = Bounded_error_squared_distance_computation<TMF, TM2, VPM1, VPM2, TMF_tree, TM2_tree, Kernel>;

    Comp bedc(tm1_parts, tm2, error_bound, vpm1, vpm2,
              infinity_value, sq_initial_bound, tm1_trees, tm2_tree);
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, tm1_parts.size()), bedc);

    sq_hdist = bedc.sq_hdist;
  }
  else // sequential version
#endif // defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_METIS_ENABLED)
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout << "* executing sequential version" << std::endl;
#endif
    sq_hdist = bounded_error_squared_Hausdorff_distance_impl<Kernel>(
                 tm1, tm2, vpm1, vpm2, tm1_tree, tm2_tree,
                 error_bound, sq_initial_bound, sq_distance_bound, infinity_value, out);
  }

#ifdef CGAL_HAUSDORFF_DEBUG
  timer.stop();
  std::cout << "* squared distance " << sq_hdist << std::endl;
  std::cout << "* distance " << approximate_sqrt(sq_hdist) << std::endl;
  std::cout << "* computation time (sec.) " << timer.time() << std::endl;
#endif

  CGAL_postcondition(sq_hdist >= FT(0));

  return sq_hdist;
}

template <class Concurrency_tag,
          class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2,
          class NamedParameters1,
          class NamedParameters2,
          class OutputIterator1,
          class OutputIterator2>
typename Kernel::FT
bounded_error_squared_symmetric_Hausdorff_distance_impl(const TriangleMesh1& tm1,
                                                        const TriangleMesh2& tm2,
                                                        const typename Kernel::FT error_bound,
                                                        const typename Kernel::FT sq_distance_bound,
                                                        const bool compare_meshes,
                                                        const VPM1 vpm1,
                                                        const VPM2 vpm2,
                                                        const NamedParameters1& np1,
                                                        const NamedParameters2& np2,
                                                        OutputIterator1& out1,
                                                        OutputIterator2& out2)
{
#if !defined(CGAL_LINKED_WITH_TBB) || !defined(CGAL_METIS_ENABLED)
  static_assert(!std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value,
                "Parallel_tag is enabled but at least TBB or METIS is unavailable.");
#endif

  // Optimized version.
  // -- We compare meshes only if it is required.
  // -- We first build trees and rebuild them only if it is required.
  // -- We provide better initial lower bound in the second call to the Hausdorff distance.
  using FT = typename Kernel::FT;

  using TM1_primitive = AABB_face_graph_triangle_primitive<TriangleMesh1, VPM1>;
  using TM2_primitive = AABB_face_graph_triangle_primitive<TriangleMesh2, VPM2>;

  using TM1_traits = AABB_traits_3<Kernel, TM1_primitive>;
  using TM2_traits = AABB_traits_3<Kernel, TM2_primitive>;

  using TM1_tree = AABB_tree<TM1_traits>;
  using TM2_tree = AABB_tree<TM2_traits>;

  using Face_handle_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using Face_handle_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;

  std::vector<Face_handle_1> tm1_only;
  std::vector<Face_handle_2> tm2_only;

  const FT sq_error_bound = square(FT(error_bound));
  FT infinity_value = FT(-1);

  // All trees below are built and/or accelerated lazily.
  TM1_tree tm1_tree;
  TM2_tree tm2_tree;
  bool rebuild = false;
  std::tie(infinity_value, rebuild) = preprocess_bounded_error_squared_Hausdorff_distance_impl<Kernel>(
    tm1, tm2, compare_meshes, vpm1, vpm2, false /*is_one_sided_distance*/, np1, np2,
    tm1_tree, tm2_tree, tm1_only, tm2_only);

  if(is_negative(infinity_value))
  {
#ifdef CGAL_HAUSDORFF_DEBUG
    std::cout.precision(17);
    std::cout << "* culling rate: 100%" << std::endl;
#endif
    const auto face1 = *(faces(tm1).begin());
    const auto face2 = *(faces(tm2).begin());
    *out1++ = std::make_pair(face1, face2);
    *out1++ = std::make_pair(face1, face2);
    *out2++ = std::make_pair(face2, face1);
    *out2++ = std::make_pair(face2, face1);

    return 0.; // TM1 and TM2 are equal so the distance is zero
  }

  CGAL_assertion(is_positive(infinity_value));

  // Compute the first one-sided distance.
  FT sq_initial_bound = sq_error_bound;
  FT sq_dista = sq_error_bound;

  if(!compare_meshes || (compare_meshes && tm1_only.size() > 0))
  {
    sq_dista = bounded_error_squared_Hausdorff_distance_impl<Kernel>(
                 tm1, tm2, vpm1, vpm2, tm1_tree, tm2_tree,
                 error_bound, sq_initial_bound, sq_distance_bound, infinity_value, out1);
  }

  // In case this is true, we need to rebuild trees in order to accelerate
  // computations for the second call.
  if(rebuild)
  {
    CGAL_assertion(compare_meshes);
    tm1_tree.clear();
    tm2_tree.clear();
    CGAL_assertion(tm2_only.size() > 0);
    CGAL_assertion(tm2_only.size() < faces(tm2).size());
    tm1_tree.insert(faces(tm1).begin(), faces(tm1).end(), tm1, vpm1);
    tm2_tree.insert(tm2_only.begin(), tm2_only.end(), tm2, vpm2);
  }

  // Compute the second one-sided distance.
  sq_initial_bound = sq_dista; // @todo we should better test this optimization!
  FT sq_distb = sq_error_bound;

  if(!compare_meshes || (compare_meshes && tm2_only.size() > 0))
  {
    sq_distb = bounded_error_squared_Hausdorff_distance_impl<Kernel>(
                 tm2, tm1, vpm2, vpm1, tm2_tree, tm1_tree,
                 error_bound, sq_initial_bound, sq_distance_bound, infinity_value, out2);
  }

  return (CGAL::max)(sq_dista, sq_distb);
}

template<class Kernel, class TM2_tree>
typename Kernel::FT recursive_hausdorff_subdivision(const typename Kernel::Point_3& p0,
                                                    const typename Kernel::Point_3& p1,
                                                    const typename Kernel::Point_3& p2,
                                                    const TM2_tree& tm2_tree,
                                                    const typename Kernel::FT sq_error_bound)
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  auto midpoint = Kernel().construct_midpoint_3_object();
  auto squared_distance = Kernel().compute_squared_distance_3_object();
  // If all edge lengths of the triangle are below the error bound,
  // return the maximum of the distances of the three points to TM2 (via TM2_tree).
  const FT max_squared_edge_length = (CGAL::max)((CGAL::max)(squared_distance(p0, p1),
                                                             squared_distance(p0, p2)),
                                                             squared_distance(p1, p2));

  if(max_squared_edge_length < sq_error_bound)
  {
    return (CGAL::max)((CGAL::max)(squared_distance(p0, tm2_tree.closest_point(p0)),
                                   squared_distance(p1, tm2_tree.closest_point(p1))),
                                   squared_distance(p2, tm2_tree.closest_point(p2)));
  }

  // Else subdivide the triangle and proceed recursively.
  const Point_3 p01 = midpoint(p0, p1);
  const Point_3 p02 = midpoint(p0, p2);
  const Point_3 p12 = midpoint(p1, p2);

  return (CGAL::max)(
           (CGAL::max)(recursive_hausdorff_subdivision<Kernel>( p0, p01, p02, tm2_tree, sq_error_bound),
                       recursive_hausdorff_subdivision<Kernel>( p1, p01, p12, tm2_tree, sq_error_bound)),
           (CGAL::max)(recursive_hausdorff_subdivision<Kernel>( p2, p02, p12, tm2_tree, sq_error_bound),
                       recursive_hausdorff_subdivision<Kernel>(p01, p02, p12, tm2_tree, sq_error_bound)));
}

template <class Concurrency_tag,
          class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2>
typename Kernel::FT
bounded_error_squared_Hausdorff_distance_naive_impl(const TriangleMesh1& tm1,
                                                    const TriangleMesh2& tm2,
                                                    const typename Kernel::FT sq_error_bound,
                                                    const VPM1 vpm1,
                                                    const VPM2 vpm2)
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  using TM2_primitive = AABB_face_graph_triangle_primitive<TriangleMesh2, VPM2>;
  using TM2_traits = AABB_traits_3<Kernel, TM2_primitive>;
  using TM2_tree = AABB_tree<TM2_traits>;

  using TM1_face_to_triangle_map = Triangle_from_face_descriptor_map<TriangleMesh1, VPM1>;

  FT sq_lower_bound = FT(0);

  // Build an AABB tree on tm2.
  TM2_tree tm2_tree(faces(tm2).begin(), faces(tm2).end(), tm2, vpm2);
  tm2_tree.build();
  tm2_tree.accelerate_distance_queries();

  // Build a map to obtain actual triangles from the face descriptors of tm1.
  const TM1_face_to_triangle_map face_to_triangle_map(&tm1, vpm1);

  // Iterate over the faces of TM1.
  for(const auto& face : faces(tm1))
  {
    // Get the vertices of the face and pass them on to a recursive method.
    const Triangle_3 triangle = get(face_to_triangle_map, face);
    const Point_3& v0 = triangle.vertex(0);
    const Point_3& v1 = triangle.vertex(1);
    const Point_3& v2 = triangle.vertex(2);

    // Recursively process the current triangle to obtain a lower bound on its Hausdorff distance.
    const FT sq_triangle_bound = recursive_hausdorff_subdivision<Kernel>(v0, v1, v2, tm2_tree, sq_error_bound);

    // Store the largest lower bound.
    if(sq_triangle_bound > sq_lower_bound)
      sq_lower_bound = sq_triangle_bound;
  }

  return to_double(approximate_sqrt(sq_lower_bound));
}

} // namespace internal

/**
 * \ingroup PMP_distance_grp
 *
 * returns an estimate on the Hausdorff distance from `tm1` to `tm2` that
 * is at most `error_bound` away from the actual Hausdorff distance from `tm1` to `tm2`.
 *
 * @tparam Concurrency_tag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag` and `Parallel_tag`.
 *                         Currently, the parallel version is not implemented and the
 *                         sequential version is always used whatever tag is chosen!
 *
 * @tparam TriangleMesh1 a model of the concept `FaceListGraph`
 * @tparam TriangleMesh2 a model of the concept `FaceListGraph`
 *
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm1 a triangle mesh
 * @param tm2 another triangle mesh
 *
 * @param error_bound a maximum bound by which the Hausdorff distance estimate is
 *                    allowed to deviate from the actual Hausdorff distance.
 *
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmX`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshX>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmX)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMeshX`.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{match_faces}
 *     \cgalParamDescription{a boolean tag that turns on the preprocessing step that filters out all faces
 *                           which belong to both meshes and hence do not contribute to the final distance}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *     \cgalParamExtra{Both `np1` and `np2` must have this tag set to `true` in order to activate this preprocessing.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @pre `tm1` and `tm2` are non-empty triangle meshes.
 *
 * @return the one-sided Hausdorff distance
 */
template <class Concurrency_tag,
          class TriangleMesh1,
          class TriangleMesh2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
double bounded_error_Hausdorff_distance(const TriangleMesh1& tm1,
                                        const TriangleMesh2& tm2,
                                        const double error_bound = 0.0001,
                                        const NamedParameters1& np1 = parameters::default_values(),
                                        const NamedParameters2& np2 = parameters::default_values())
{
  using Traits = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  using FT = typename Traits::FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(!is_empty(tm1) && is_triangle_mesh(tm1));
  CGAL_precondition(!is_empty(tm2) && is_triangle_mesh(tm2));

  const auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm1));
  const auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm2));

  const bool match_faces1 = choose_parameter(get_parameter(np1, internal_np::match_faces), true);
  const bool match_faces2 = choose_parameter(get_parameter(np2, internal_np::match_faces), true);
  const bool match_faces = match_faces1 && match_faces2;

  auto out = choose_parameter(get_parameter(np1, internal_np::output_iterator),
                              CGAL::Emptyset_iterator());

  CGAL_precondition(error_bound >= 0.);

  const FT sq_hdist = internal::bounded_error_squared_one_sided_Hausdorff_distance_impl<Concurrency_tag, Traits>(
        tm1, tm2, error_bound, FT(-1) /*distance threshold*/, match_faces, vpm1, vpm2, np1, np2, out);

  return to_double(approximate_sqrt(sq_hdist));
}

/**
 * \ingroup PMP_distance_grp
 *
 * returns the the symmetric Hausdorff distance, that is
 * the maximum of `bounded_error_Hausdorff_distance(tm1, tm2, error_bound, np1, np2)`
 * and `bounded_error_Hausdorff_distance(tm2, tm1, error_bound, np2, np1)`.
 *
 * This function optimizes all internal calls to shared data structures in order to
 * speed up the computation.
 *
 * See the function `CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance()`
 * for a complete description of the parameters and requirements.
 */
template <class Concurrency_tag,
          class TriangleMesh1,
          class TriangleMesh2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
double bounded_error_symmetric_Hausdorff_distance(const TriangleMesh1& tm1,
                                                  const TriangleMesh2& tm2,
                                                  const double error_bound,
                                                  const NamedParameters1& np1 = parameters::default_values(),
                                                  const NamedParameters2& np2 = parameters::default_values())
{
  using Traits = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  using FT = typename Traits::FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(!is_empty(tm1) && is_triangle_mesh(tm1));
  CGAL_precondition(!is_empty(tm2) && is_triangle_mesh(tm2));

  const auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm1));
  const auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm2));

  const bool match_faces1 = choose_parameter(get_parameter(np1, internal_np::match_faces), true);
  const bool match_faces2 = choose_parameter(get_parameter(np2, internal_np::match_faces), true);
  const bool match_faces = match_faces1 && match_faces2;

  // TODO: should we return a union of these realizing triangles?
  auto out1 = choose_parameter(get_parameter(np1, internal_np::output_iterator),
                               CGAL::Emptyset_iterator());
  auto out2 = choose_parameter(get_parameter(np2, internal_np::output_iterator),
                               CGAL::Emptyset_iterator());

  CGAL_precondition(error_bound >= 0.);

  const FT sq_hdist = internal::bounded_error_squared_symmetric_Hausdorff_distance_impl<Concurrency_tag, Traits>(
    tm1, tm2, error_bound, FT(-1) /*distance_threshold*/, match_faces, vpm1, vpm2, np1, np2, out1, out2);

  return to_double(approximate_sqrt(sq_hdist));
}

/**
 * \ingroup PMP_distance_grp
 *
 * \brief returns `true` if the Hausdorff distance between two meshes is larger than
 * the user-defined max distance, otherwise it returns `false`.
 *
 * The distance used to compute the proximity of the meshes is the bounded-error Hausdorff distance.
 * Instead of computing the full distance and checking it against the user-provided
 * value, this function returns early if certain criteria show that the meshes
 * do not satisfy the provided `distance_bound`.
 *
 * See the function `CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance()`
 * for a complete description of the parameters and requirements. The following extra named parameter
 * is available for `np1`:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_one_sided_hausdorff}
 *     \cgalParamDescription{a boolean tag indicating if the one-sided Hausdorff distance should be used.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *     \cgalParamExtra{If this tag is set to `false`, the symmetric Hausdorff distance is used.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template< class Concurrency_tag,
          class TriangleMesh1,
          class TriangleMesh2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
bool is_Hausdorff_distance_larger(const TriangleMesh1& tm1,
                                  const TriangleMesh2& tm2,
                                  const double distance_bound,
                                  const double error_bound,
                                  const NamedParameters1& np1 = parameters::default_values(),
                                  const NamedParameters2& np2 = parameters::default_values())
{
  using Traits = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  using FT = typename Traits::FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(!is_empty(tm1) && is_triangle_mesh(tm1));
  CGAL_precondition(!is_empty(tm2) && is_triangle_mesh(tm2));

  if(distance_bound <= 0.)
    return true;

  const auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm1));
  const auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm2));

  const bool match_faces1 = choose_parameter(get_parameter(np1, internal_np::match_faces), true);
  const bool match_faces2 = choose_parameter(get_parameter(np2, internal_np::match_faces), true);
  const bool match_faces = match_faces1 && match_faces2;
  const bool use_one_sided = choose_parameter(get_parameter(np1, internal_np::use_one_sided_hausdorff), true);

  CGAL_precondition(error_bound >= 0.);
  CGAL_precondition(distance_bound > 0.);

  const FT sq_distance_bound = square(FT(distance_bound));

  auto stub = CGAL::Emptyset_iterator();

  FT sq_hdist = FT(-1);
  if(use_one_sided)
  {
    sq_hdist = internal::bounded_error_squared_one_sided_Hausdorff_distance_impl<Concurrency_tag, Traits>(
      tm1, tm2, error_bound, sq_distance_bound, match_faces, vpm1, vpm2, np1, np2, stub);
  }
  else
  {
    sq_hdist = internal::bounded_error_squared_symmetric_Hausdorff_distance_impl<Concurrency_tag, Traits>(
      tm1, tm2, error_bound, sq_distance_bound, match_faces, vpm1, vpm2, np1, np2, stub, stub);
  }

#ifdef CGAL_HAUSDORFF_DEBUG
  std::cout.precision(17);
  std::cout << "- fin distance: " << approximate_sqrt(sq_hdist) << std::endl;
  std::cout << "- max distance: " << distance_bound << std::endl;
#endif

  return (sq_hdist > sq_distance_bound);
}

// Implementation of the naive Bounded Error Hausdorff distance.
template <class Concurrency_tag,
          class TriangleMesh1,
          class TriangleMesh2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
double bounded_error_Hausdorff_distance_naive(const TriangleMesh1& tm1,
                                              const TriangleMesh2& tm2,
                                              const double error_bound,
                                              const NamedParameters1& np1 = parameters::default_values(),
                                              const NamedParameters2& np2 = parameters::default_values())
{
  using Traits = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  using FT = typename Traits::FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(!is_empty(tm1) && is_triangle_mesh(tm1));
  CGAL_precondition(!is_empty(tm2) && is_triangle_mesh(tm2));

  const auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm1));
  const auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                     get_const_property_map(vertex_point, tm2));

  CGAL_precondition(error_bound >= 0.);

  const FT sq_hdist = internal::bounded_error_squared_Hausdorff_distance_naive_impl<Concurrency_tag, Traits>(
           tm1, tm2, error_bound, vpm1, vpm2);

  return to_double(approximate_sqrt(sq_hdist));
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
